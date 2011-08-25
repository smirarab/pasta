#!/usr/bin/env python

"""Multi-threaded jobs
"""

# This file is part of SATe

# SATe is free software: you can redistribute it and/or modify
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
from cStringIO import StringIO
from Queue import Queue
from threading import Thread, Event, Lock
from subprocess import Popen, PIPE
from sate import get_logger, TIMING_LOG
import random, sys
from sate.filemgr import open_with_intermediates

_LOG = get_logger(__name__)

class LoggingQueue(Queue):
    def put(self, job):
        TIMING_LOG.info("%s queued" % str(job.context_str))
        Queue.put(self, job)

jobq = LoggingQueue()

def worker():
    while True:
        job = jobq.get()
        TIMING_LOG.info("%s started" % str(job.context_str))
        try:
            job.start()
        except:
            err = StringIO()
            traceback.print_exc(file=err)
            _LOG.error("Worker dying.  Error in job.start = %s" % err.getvalue())
        else:
            try:
                job.get_results()
            except:
                err = StringIO()
                traceback.print_exc(file=err)
                _LOG.error("Worker dying.  Error in job.get_results = %s" % err.getvalue())
        TIMING_LOG.info("%s completed" % str(job.context_str))
        jobq.task_done()
    return

# We'll keep a list of Worker threads that are running in case any of our code triggers multiple calls
_WORKER_THREADS = []

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
        t = Thread(target=worker)
        _WORKER_THREADS.append(t)
        t.setDaemon(True)
        t.start()

class JobBase(object):
    def __init__(self, **kwargs):
        self.context_str = kwargs.get("context_str")
        if "context_str" in kwargs:
            del kwargs['context_str']
        self._kwargs = kwargs

class FakeJob(JobBase):
    """FakeJob instances are used in cases in which we know have the results of
    an operation, but need to emulate the API of DispatchableJob.
    """
    def __init__(self, results, **kwargs):
        JobBase.__init__(self, **kwargs)
        self.results = results

    def start(self):
        pass

    def wait(self):
        pass

    def get_results(self):
        return self.results

    def kill(self):
        pass


class DispatchableJob(JobBase):
    def __init__(self, invocation, result_processor, **kwargs):
        JobBase.__init__(self, **kwargs)
        self._invocation = invocation
        # _LOG.debug('DispatchableJob.__init__(invocation= %s )' % " ".join(self._invocation))  # Not sure why it does not work with datatype in treebuild.create_job
        self.result_processor = result_processor
        self.process = None
        self.return_code = None
        self.results = None
        self._id = None
        self._stdout_fo = None
        self._stderr_fo = None
        self.launched_event = Event()
        self.finished_event = Event()
        self.thread_waiting = False
        self.thread_waiting_lock = Lock()
        self.error = None

    def get_id(self):
        return self._id

    def set_id(self, i):
        self._id = i

    id = property(get_id, set_id)

    def start(self):
        assert self.process is None, "Relaunching jobs is not allowed"
        try:
            _LOG.debug('launching %s.\n setting event' % " ".join(self._invocation))
            proc_cwd = self._kwargs.get('cwd', os.curdir)
            k = dict(self._kwargs)
            if 'stdout' not in self._kwargs:
                self._stdout_fo = open_with_intermediates(os.path.join(proc_cwd, '.Job.stdout.txt'), 'w')
                k['stdout'] = self._stdout_fo
            if 'stderr' not in self._kwargs:
                self._stderr_fo = open_with_intermediates(os.path.join(proc_cwd, '.Job.stderr.txt'), 'w')
                k['stderr'] = self._stderr_fo
            self.process = Popen(self._invocation, stdin = PIPE, **k)
            self.set_id(self.process.pid)
            #f = open('.%s.pid' % self.get_id(), 'w')
            #f.close()
            _LOG.debug('setting launched_event')
        except:
            self.error = RuntimeError('The invocation:\n"%s"\nfailed' % '" "'.join(self._invocation))
            raise
        finally:
            self.launched_event.set()

####
#   Polling does not appear to be needed in the current impl.
#   def poll(self):
#       "Not blocking"
#       if self.return_code is None:
#           self.return_code = self.process.poll()
#       return self.return_code

    def read_stderr(self):
        if os.path.exists(self._stderr_fo.name):
            errFile = open(self._stderr_fo.name,'r')
            errorFromFile = errFile.read(-1)
            errFile.close()
            return errorFromFile
        else:
            return None

    def wait(self):
        """Blocking.

        I'm not sure that it is safe for multiple threads calling self.process.wait
        so we'll only have the first thread do this.  All other threads that enter
        wait will wait for the finished_event
        """
        self.thread_waiting_lock.acquire()
        if self.thread_waiting:
            self.thread_waiting_lock.release()

            if self.error is not None:
                raise self.error


            self.finished_event.wait()
            if self.error is not None:
                raise self.error
        else:
            # this branch is actually monitoring the process
            self.thread_waiting = True
            self.thread_waiting_lock.release()

            try:
                if self.error is not None:
                    raise self.error

                if self.return_code is None:
                    _LOG.debug('Launch detected')
                    self.launched_event.wait()
                    _LOG.debug('Launch detected')
                    if self.error is not None:
                        raise self.error

                    err_msg = []
                    err_msg.append("SATe failed because one of the programs it tried to run failed.")
                    err_msg.append('The invocation that failed was: \n    "%s"\n' % '" "'.join(self._invocation))

                    try:
                        self.return_code = self.process.wait()
                        if self._stdout_fo:
                            self._stdout_fo.close()
                        if self._stderr_fo:
                            self._stderr_fo.close()
                        if self.return_code:
                            errorFromFile = self.read_stderr()
                            if errorFromFile:
                                err_msg.append(errorFromFile)
                            self.error = Exception("\n".join(err_msg))
                            raise self.error
                    except:
                        self.error = Exception("\n".join(err_msg))
                        raise self.error

                    try:
                        self.results = self.result_processor()
                    except Exception, e:
                        errorFromFile = self.read_stderr()
                        if not errorFromFile:
                            errorFromFile = str(e)
                        err_msg.append(errorFromFile)
                        self.error = Exception("\n".join(err_msg))

                    if self.error is not None:
                        raise self.error
            finally:
                self.finished_event.set()

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

    def kill(self):
        if self.results is None:
            self.process.kill()
