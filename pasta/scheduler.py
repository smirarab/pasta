#!/usr/bin/env python

"""Multi-threaded jobs
"""

# This file is part of PASTA and is forked from SATe

# PASTA like SATe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Jiaye Yu and Mark Holder, University of Kansas

import os, traceback
from io import StringIO
from io import BytesIO
try:
    from queue import Queue
except ImportError:
    from Queue import Queue
from threading import Thread, Event, Lock
from multiprocessing import Process, Manager, Value
from subprocess import Popen, PIPE
from pasta import get_logger, TIMING_LOG
from pasta.filemgr import open_with_intermediates
from random import random

_LOG = get_logger(__name__)

class LoggingQueue(Queue):
    def put(self, job):
        TIMING_LOG.info("%s queued" % str(job.context_str))
        _LOG.debug("%s queued" % str(job.context_str))
        Queue.put(self, job)

jobq = LoggingQueue()

_all_dispatchable_jobs = []

merged_queue_events = []

def new_merge_event():
    global merged_queue_events
    e = Event()
    merged_queue_events.append(e)
    return e

def set_all_events():
    global merged_queue_events
    for e in merged_queue_events:
        e.set()
        
def kill_all_jobs():
    global _all_dispatchable_jobs
    for job in _all_dispatchable_jobs:
        job.kill()
    set_all_events()

class LightJobForProcess():

    def __init__(self, invocation, k, environ):
        self._invocation = invocation
        self._k = k
        self.error = None # NOTE: This is NOT sent back to the main process
        self.return_code = None # NOTE: This is sent back to the main process
        self.environ = environ

    def read_stderr(self,_stderr_fo):
        if os.path.exists(_stderr_fo.name):
            errFile = open(_stderr_fo.name,'r')
            errorFromFile = errFile.read(-1)
            errFile.close()
            return errorFromFile
        else:
            return None
                
    def run(self):
        _LOG.debug('launching %s.' % " ".join(self._invocation))
        k = self._k
        proc_cwd = k.get('cwd', os.curdir)
        stdout_file_path = k.get('stdout', None)
        stderr_file_path = k.get('stderr', None)
        if stdout_file_path:
            _stdout_fo = open_with_intermediates(stdout_file_path, 'w')
        else:
            _stdout_fo = open_with_intermediates(os.path.join(proc_cwd, '.Job.stdout.txt'), 'w')
        k['stdout'] = _stdout_fo
        if stderr_file_path:
            _stderr_fo = open_with_intermediates(stderr_file_path, 'w')
        else:
            _stderr_fo = open_with_intermediates(os.path.join(proc_cwd, '.Job.stderr.txt'), 'w')
        k['stderr'] = _stderr_fo

        for key,v in self.environ.items():
            os.environ[key] = v

        process = Popen(self._invocation, stdin = PIPE, universal_newlines=True, **k)

        err_msg = []                
        err_msg.append("PASTA failed because one of the programs it tried to run failed.")
        err_msg.append('The invocation that failed was: \n    "%s"\n' % '" "'.join(self._invocation))
        try:
            self.return_code = process.wait()
            _stdout_fo.close()
            _stderr_fo.close()
            process.stdin.close()
            if self.return_code:
                errorFromFile = self.read_stderr(_stderr_fo)
                if errorFromFile:
                    err_msg.append(errorFromFile)
                self.error = "\n".join(err_msg)
                raise Exception("")
            _LOG.debug('Finished %s.\n Return code: %s; %s' % (" ".join(self._invocation), self.return_code, self.error))        
        except Exception as e:
            err_msg.append(str(e))
            self.error = "\n".join(err_msg) 
            _LOG.error(self.error)   
            
class pworker():
    def __init__(self, i, q, err_shared_obj):
        self.q = q
        self.i = i
        self.err_shared_obj = err_shared_obj
        
    def __call__(self):
        while True:                                
            try:
                _LOG.debug("Process Worker %d ticking on queue %s of size %d" %(self.i,str(self.q),self.q.qsize()))
                job = self.q.get()
                _LOG.debug("Process Worker %d found a job to run" %self.i)
                plj = LightJobForProcess(job[0],job[1],job[2])
                plj.run()
                self.err_shared_obj.value = plj.return_code
                job[3] = plj.error
                self.q.task_done()
            except:
                err = BytesIO()
                traceback.print_exc(file=err)
                _LOG.error("Process Worker dying.  Error in job.start = %s" % err.getvalue())
                raise
        return

_manager = Manager()

class worker():
    
    def __init__(self, i):
        global _manager
        self.i = i
        self.pqueue = _manager.Queue()
        self.err_shared_obj = Value('i', 0)
        pw = pworker(self.i, self.pqueue, self.err_shared_obj)
        self.p = Process(target=pw)
        self.p.daemon = True
        self.p.start()

    def stop(self):
        self.p.terminate()
        
    def __call__(self):                            
        while True:            
            job = jobq.get()            
            ID=int(random() *10000000)
            TIMING_LOG.info("%s (%d) started" % (str(job.context_str),ID))
            try:
                if isinstance(job, DispatchableJob):
                    pa = job.start()
                    shared_job_obj = [pa[0],pa[1],dict(os.environ),None]
                    self.pqueue.put(shared_job_obj)
                    _LOG.debug("Worker %d put a job tuple on queue %s" %(self.i,str(self.pqueue)))
                    
                    self.pqueue.join()   

                    _LOG.debug("Worker %d joined on queue %s" %(self.i,str(self.pqueue)))
                    
                    plj = LightJobForProcess(shared_job_obj[0],shared_job_obj[1],shared_job_obj[2])
                    plj.error = shared_job_obj[3]
                    plj.return_code = self.err_shared_obj.value
                                             
                    if plj.error is not None:
                        job.error = Exception(plj.error)
                        
                    job.return_code = plj.return_code
                    
                    if job.return_code is not None and job.return_code != 0:
                        raise Exception("Job:\n %s\n failed with error code: %d" %(' '.join(plj._invocation), job.return_code))
                    
                    job.results = job.result_processor()
                    
                    job.finished_event.set() 
                    job.get_results()
                    job.postprocess()
                else:                    
                    job.start()
                    job.get_results()
                    job.postprocess()

            except Exception as e:
                err = BytesIO()
                traceback.print_exc(file=err)
                _LOG.error("Worker dying.  Error in job.start = %s" % err.getvalue())
                job.error=e
                job.return_code = -1
                job.finished_event.set() 
                job.kill()
                kill_all_jobs()
                return                
            TIMING_LOG.info("%s (%d) completed" % (str(job.context_str),ID))
            jobq.task_done()
        return

# We'll keep a list of Worker threads that are running in case any of our code triggers multiple calls
_WORKER_THREADS = []
_WORKER_OBJECTS = []

def start_worker(num_workers):
    """Spawns worker threads such that at least `num_workers` threads will be
    launched for processing jobs in the jobq.

    The only way that you can get more than `num_workers` threads is if you
    have previously called the function with a number > `num_workers`.
    (worker threads are never killed).
    """
    assert num_workers > 0, "A positive number must be passed as the number of worker threads"
    num_currently_running = len(_WORKER_THREADS)
    for i in range(num_currently_running, num_workers):
        _LOG.debug("Launching Worker thread #%d" % i)
        w = worker(i)
        t = Thread(target=w)
        _WORKER_THREADS.append(t)
        _WORKER_OBJECTS.append(w)
        t.setDaemon(True)
        t.start()
        
def stop_worker():
    for w in _WORKER_OBJECTS:
        w.stop()

class JobBase(object):
    def __init__(self, **kwargs):
        self.context_str = kwargs.get("context_str")
        if "context_str" in kwargs:
            del kwargs['context_str']
        self._kwargs = kwargs


class DispatchableJob(JobBase):
    def __init__(self, invocation, result_processor, **kwargs):
        global _all_dispatchable_jobs
        JobBase.__init__(self, **kwargs)
        self._invocation = invocation
        # _LOG.debug('DispatchableJob.__init__(invocation= %s )' % " ".join(self._invocation))  # Not sure why it does not work with datatype in treebuild.create_job
        self.result_processor = result_processor
        self.return_code = None
        self.results = None
        self._id = None
        self._stdout_fo = None
        self._stderr_fo = None
        self.finished_event = Event()        
        self.error = None
        _all_dispatchable_jobs.append(self)

    def get_id(self):
        return self._id

    def set_id(self, i):
        self._id = i

    id = property(get_id, set_id)

    def start(self):
        try:
            #_LOG.debug('launching %s.\n setting event' % " ".join(self._invocation))            
            k = dict(self._kwargs)
            return (self._invocation, k)
            #self.process = Popen(self._invocation, stdin = PIPE, **k)
            #self.set_id(self.process.pid)
            #f = open('.%s.pid' % self.get_id(), 'w')
            #f.close()            
        except:
            self.error = RuntimeError('The invocation:\n"%s"\nfailed' % '" "'.join(self._invocation))
            raise

####
#   Polling does not appear to be needed in the current impl.
#   def poll(self):
#       "Not blocking"
#       if self.return_code is None:
#           self.return_code = self.process.poll()
#       return self.return_code

    def wait(self):
        """Blocking.

        I'm not sure that it is safe for multiple threads calling self.process.wait
        so we'll only have the first thread do this.  All other threads that enter
        wait will wait for the finished_event
        """
        
        # this branch is actually monitoring the process
        if self.error is not None:
            raise self.error

        self.finished_event.wait() 
       
        if self.error is not None:
            raise self.error                        
        
        return self.return_code

    def get_results(self):
        if self.error is not None:
            raise self.error
        if self.results is None:
            self.wait()
        else:
            pass
            # os.remove('.%s.pid' % self.get_id())  # treebuild_job will not go through this after finish
        return self.results
    
    def postprocess(self):
        pass
    
    def kill(self):
        self.finished_event.set()
            
            
class TickableJob():

    def __init__(self):
        self._parents = []
        self._unfinished_children = []
        self._childrenlock = Lock()
        self._parentslock = Lock()

    def on_dependency_ready(self):
        ''' This function needs to be implemented by the children. 
                    This will be called when all children are done'''
        raise NotImplementedError("on_dependency_ready() method is not implemented.")

    def add_parent(self,parentjob):
        self._parentslock.acquire()
        self._parents.append(parentjob)
        self._parentslock.release()	

    def add_child(self,childJob):
        self._childrenlock.acquire()
        self._unfinished_children.append(childJob)
        self._childrenlock.release()

    def tick(self, finishedjob):
        #_LOG.debug("ticking %s" %str(self))
        self._childrenlock.acquire()
        self._unfinished_children.remove(finishedjob)   
        #_LOG.debug("children ... %d" %len(self._unfinished_children))     
        nochildleftbehind = not bool(self._unfinished_children)
        self._childrenlock.release()
        if nochildleftbehind:
            self.on_dependency_ready()
            
    def tick_praents(self):
        self._parentslock.acquire()
        for parent in self._parents:
            parent.tick(self)
        self._parentslock.release()  
        
    def kill(self):    
        self._parentslock.acquire()
        for parent in self._parents:
            parent.kill()
        self._parentslock.release()          

class TickingJob():
    def __init__(self):
        self.parent_tickable_job = []
        self.killed = False
    
    def add_parent_tickable_job(self, tickableJob):
        self.parent_tickable_job.append(tickableJob)
        
    def postprocess(self):
        for parent in self.parent_tickable_job:
            parent.tick(self)
            
    def kill(self):
        if not self.killed:
            self.killed = True
            for parent in self.parent_tickable_job:
                parent.kill()            

class TickingDispatchableJob(DispatchableJob,TickingJob):
    def __init__(self, invocation, result_processor, **kwargs):
        DispatchableJob.__init__(self, invocation, result_processor, **kwargs)
        TickingJob.__init__(self)

    def kill(self):
        DispatchableJob.kill(self)
        TickingJob.kill(self)

    def postprocess(self):
        TickingJob.postprocess(self)


class FakeJob(JobBase, TickingJob):
    """FakeJob instances are used in cases in which we know have the results of
    an operation, but need to emulate the API of DispatchableJob.
    """
    def __init__(self, results, **kwargs):
        JobBase.__init__(self, **kwargs)
        TickingJob.__init__(self)
        self.results = results

    def start(self):
        pass

    def wait(self):
        pass

    def get_results(self):
        return self.results
    
    def kill(self):
        pass
