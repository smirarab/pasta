"""Main script for PASTA GUI on Windows/Mac/Linux
"""

# This file is part of PASTA which is forked from SATe

# PASTA, like SATe is free software: you can redistribute it and/or modify
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

import os
import platform
import subprocess
import tempfile
import sys
import time
import wx
import string
from pasta import PROGRAM_AUTHOR
from pasta import PROGRAM_INSTITUTE
from pasta import PROGRAM_DESCRIPTION
from pasta import PROGRAM_LICENSE
from pasta import PROGRAM_NAME
from pasta import PROGRAM_VERSION
from pasta import PROGRAM_WEBSITE
from pasta import PROGRAM_YEAR
from pasta import GLOBAL_DEBUG
from pasta import DEFAULT_MAX_MB
try:
    from configparser import RawConfigParser
except:
    from ConfigParser import RawConfigParser
from pasta import pasta_is_frozen
from pasta import pasta_home_dir
from pasta.configure import get_invoke_run_pasta_command
from pasta.tools import AlignerClasses
from pasta.tools import MergerClasses
from pasta.tools import TreeEstimatorClasses
from pasta.tools import get_aligner_classes, get_merger_classes, get_tree_estimator_classes
from pasta import filemgr
from pasta.usersettingclasses import get_list_of_seq_filepaths_from_dir
from pasta.alignment import summary_stats_from_parse
from pasta.mainpasta import get_auto_defaults_from_summary_stats

WELCOME_MESSAGE = "%s %s, %s\n\n"% (PROGRAM_NAME, PROGRAM_VERSION, PROGRAM_YEAR)
GRID_VGAP = 8
GRID_HGAP = 8

PARSING_FILES_IN_GUI = True
MAX_NUM_CPU = 16
PASTA_GUI_ONLY_PRINTS_CONFIG = os.environ.get('PASTA_GUI_ONLY_PRINTS_CONFIG') == '1'

def is_valid_int_str(s, min_v, max_v):
    try:
        i = int(s)
    except:
        return False
    si = str(i)
    if si != s:
        return False
    if min_v is not None and i < min_v:
        return False
    if max_v is not None and i > max_v:
        return False
    return True

class RangedIntValidator(wx.PyValidator):
    def __init__(self, min_v, max_v):
        wx.PyValidator.__init__(self)
        self.min_v = min_v
        self.max_v = max_v
        self.Bind(wx.EVT_CHAR, self.OnChar)
    def Clone(self):
        return RangedIntValidator(self.min_v, self.max_v)
    def is_valid_str(self, s):
        return is_valid_int_str(s, self.min_v, self.max_v)
    def Validate(self, win):
        v = win.GetValue()
        return self.is_valid_str(v)
    def OnChar(self, event):
        key = event.GetKeyCode()
        textCtrl = self.GetWindow()
        if key == wx.WXK_BACK or key == wx.WXK_DELETE:
            textCtrl.SetBackgroundColour("white")
            event.Skip()
            return
        if key < wx.WXK_SPACE or key > 255:
            textCtrl.SetBackgroundColour("white")
            event.Skip()
            return
        if chr(key) in string.digits:    
            textCtrl.SetBackgroundColour("white")
            event.Skip()
            return
        if not wx.Validator_IsSilent():
            wx.Bell()
        # Returning without calling even.Skip eats the event before it
        # gets to the text control
        return
    def TransferToWindow(self):
         return True
    def TransferFromWindow(self):
         return True

class PastaFrame(wx.Frame):
    def __init__(self, size):
        wx.Frame.__init__(self, None, -1, "PASTA - Practical Alignment using SATe and TraAnsitivity", size=(640,480), style=wx.DEFAULT_FRAME_STYLE)
        self.SetBackgroundColour(wx.LIGHT_GREY)
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText("PASTA Ready!")
        if wx.Platform == "__WXMSW__" or wx.Platform == "__WXMAC__":
            import base64
            import io
            icon = wx.EmptyIcon()
            icon.CopyFromBitmap(wx.BitmapFromImage(wx.ImageFromStream(io.BytesIO(base64.b64decode(ICO_STR)))))
            self.SetIcon(icon)

        self.ctrls = []
        sizer_all = wx.BoxSizer(wx.VERTICAL)

        self.sizer_tool_settings = self._create_tools_sizer()
        self.sizer_data = self._create_data_sizer()
        self.sizer_pasta_settings = self._create_pasta_settings_sizer()
        self.sizer_job_settings = self._create_job_settings_sizer()
        self.sizer_workflow_settings = self._create_workflow_settings_sizer()

        sizer1 = wx.BoxSizer(wx.VERTICAL)
        sizer1.Add(self.sizer_tool_settings, 0, wx.EXPAND|wx.BOTTOM|wx.RIGHT, 5)
        sizer1.Add(self.sizer_data, 0, wx.EXPAND|wx.TOP|wx.RIGHT, 5)
        sizer1.Add(self.sizer_workflow_settings, 0, wx.EXPAND|wx.TOP|wx.RIGHT, 5)
        self.sizer_settings = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_settings.Add(sizer1, 0, wx.EXPAND|wx.ALL, 0)

        sizer2 = wx.BoxSizer(wx.VERTICAL)
        sizer2.Add(self.sizer_job_settings, 0, wx.EXPAND|wx.ALL, 0)
        sizer2.Add(self.sizer_pasta_settings, 0, wx.EXPAND|wx.ALL, 0)
        self.sizer_settings.Add(sizer2, 0, wx.EXPAND|wx.ALL, 0)

        sizer_all.Add(self.sizer_settings, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 10)

        self.button = wx.Button(self, label="Start")
        self.log = wx.TextCtrl(self, -1, "", size=(200,120),style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH2)
        self.log.AppendText(WELCOME_MESSAGE)
        self.log.AppendText("Running Log (%s %s)\n\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), time.tzname[0]))
        sizer_all.Add(self.button, 0, wx.BOTTOM|wx.ALIGN_CENTER, 10)
        sizer_all.Add(self.log, 4, wx.EXPAND)

        self.SetAutoLayout(True)
        self.Layout()
        self.SetSizerAndFit(sizer_all)

        self._create_menu()
        self.process = None
        self.process_cfg_file = None

        self.Bind(wx.EVT_IDLE, self.OnIdle)
        self.Bind(wx.EVT_END_PROCESS, self.OnProcessEnded)
        self.Bind(wx.EVT_BUTTON, self.OnButton, self.button)
        
        self.set_char_model() # this fixes the model based on the current default tree estimator


    def _create_job_settings_sizer(self):
        staticboxsizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Job Settings"), wx.VERTICAL)
        sizer = wx.GridBagSizer(GRID_VGAP, GRID_HGAP)
        cr = 0
        sizer.Add(wx.StaticText(self, -1, "Job Name"),(cr,0), flag=wx.ALIGN_LEFT )
        self.txt_jobname = wx.TextCtrl(self,-1,"pastajob")
        sizer.Add(self.txt_jobname, (cr,1), flag=wx.EXPAND)
        cr += 1
        self.outputdir_btn = wx.Button(self, label="Output Dir." )
        sizer.Add(self.outputdir_btn,(cr,0), flag=wx.ALIGN_LEFT )
        self.txt_outputdir = wx.TextCtrl(self, -1, "", size=(250,9))
        sizer.Add(self.txt_outputdir, (cr,1), flag=wx.EXPAND)
        cr += 1
        sizer.Add(wx.StaticText(self, -1, "CPU(s) Available"), (cr,0), flag=wx.ALIGN_LEFT )
        self.cb_ncpu = wx.ComboBox(self, -1, "1", choices=list(map(str, list(range(1, MAX_NUM_CPU + 1)))), style=wx.CB_READONLY)
        sizer.Add(self.cb_ncpu, (cr,1), flag=wx.EXPAND)
        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Max. Memory (MB)"), (cr,0), flag=wx.ALIGN_LEFT )
        self.txt_maxmb = wx.TextCtrl(self, -1, str(DEFAULT_MAX_MB))
        sizer.Add(self.txt_maxmb, (cr,1), flag=wx.EXPAND)

        staticboxsizer.Add(sizer, 0, wx.CENTER, 0)
        self.Bind(wx.EVT_BUTTON, self.OnChooseOutputDir, self.outputdir_btn)
        return staticboxsizer

    def validate_max_mb(self, value):
        try:
            mb = int(value)
            if mb <= 0:
                raise ValueError
            return True
        except ValueError:
            wx.MessageBox("Invalid value for Maximum MB: '" + value + "': require positive integer value.",
                    "Invalid Value for Maximum MB",
                    wx.OK|wx.ICON_EXCLAMATION)
            return False

    def OnChooseOutputDir(self, event):
        dialog = wx.DirDialog(None, "Choose directory for output", style=wx.FD_OPEN)
        dialog.ShowModal()
        self.txt_outputdir.SetValue( dialog.GetPath() )

    def _set_custom_pasta_settings(self, event):
        #self.cb_sate_presets.SetValue("(custom)")
        pass

    def _create_tools_sizer(self):
        from pasta.configure import get_configuration
        cfg = get_configuration()
        staticboxsizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "External Tools"), wx.VERTICAL)
        sizer = wx.FlexGridSizer(0, 2, GRID_VGAP, GRID_HGAP)
        items = ["Aligner", "Merger", "TreeEstimator"]
        tool_list_list = [get_aligner_classes(), get_merger_classes(), get_tree_estimator_classes()]
        self.raxml_dna_models = ["GTRCAT", "GTRGAMMA", "GTRGAMMAI"]
        self.fasttree_dna_models = ["GTR+G20", "GTR+CAT", "JC+G20", "JC+CAT"]
        prot_matrix = ["DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG"]
        prot_type = ["PROTCAT", "PROTCATI", "PROTGAMMA", "PROTGAMMAI"]

        self.raxml_prot_models = [j+i for i in prot_matrix for j in prot_type]
        self.raxml_prot_models.extend([j+i+"F" for i in prot_matrix for j in prot_type])
        self.fasttree_prot_models = ["JTT+G20", "JTT+CAT", "WAG+G20", "WAG+CAT"]

        if GLOBAL_DEBUG:
            defaults = {"Aligner":"PADALIGNER", "Merger":"PADALIGNER", "TreeEstimator":"RANDTREE"}
        else:
            defaults = {"Aligner":"MAFFT", "Merger":"MUSCLE", "TreeEstimator":"FASTTREE"}
        self.cb_tools = {}
        for item_idx, item in enumerate(items):
            text = wx.StaticText(self, -1, "Tree Estimator") if item == "TreeEstimator" else wx.StaticText(self, -1, item)
            sizer.Add(text, 0, wx.LEFT)
            tool_list = tool_list_list[item_idx]
            active_tool_name_list = []
            for tool in tool_list:
                try:
                    tool_attr_name = tool.section_name.split()[0].lower()
                    tool_path = getattr(cfg, tool_attr_name).path
                    if os.path.exists(tool_path):
                        active_tool_name_list.append(tool_attr_name.upper())
                    else:
                        print("WARNING: tool %s not found at %s" %(tool,tool_path))
                except :
                    raise
            combobox = wx.ComboBox(self, -1, defaults[item], (-1,-1), (-1,-1), active_tool_name_list, wx.CB_READONLY)
            self.cb_tools[item.lower()] = combobox
            self.ctrls.append(self.cb_tools[item.lower()])
            sizer.Add(combobox, 0, wx.EXPAND)

        self.Bind(wx.EVT_COMBOBOX, self.OnTreeEstimatorChange, self.cb_tools["treeestimator"])

        combobox = wx.ComboBox(self, -1, "GTRCAT", (-1,-1), (-1,-1), self.raxml_dna_models, wx.CB_READONLY)
        self.cb_tools["model"] = combobox
        self.ctrls.append(self.cb_tools["model"])
        sizer.Add(wx.StaticText(self, -1, "Model"), wx.LEFT)
        sizer.Add(combobox, 0, wx.EXPAND)
        staticboxsizer.Add(sizer, 0, wx.CENTER, 0)
        return staticboxsizer

    def _create_data_sizer(self):
        staticboxsizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Sequences and Tree"), wx.VERTICAL)
        sizer = wx.FlexGridSizer(0, 2, GRID_VGAP, GRID_HGAP)
        self.datatype = wx.ComboBox(self, -1, "DNA", (-1, -1), (-1, -1), ["DNA", "RNA", "Protein"], wx.CB_READONLY)
        self.seq_btn = wx.Button(self, label="Sequence file ..." )
        self.tree_btn = wx.Button(self, label="Tree file (optional) ..." )

        self.txt_seqfn = wx.TextCtrl(self,-1)
        self.txt_treefn = wx.TextCtrl(self,-1)

        self.cb_multilocus = wx.CheckBox(self, -1, "Multi-Locus Data")
        self.cb_multilocus.Disable()
        self.checkbox_aligned = wx.CheckBox(self, -1, "Use for inital tree")
        self.checkbox_aligned.SetValue(False)
        self._could_be_aligned = False
        self.checkbox_aligned.Disable()

        sizer.AddMany([ (self.seq_btn, 0, wx.LEFT|wx.EXPAND),
                        (self.txt_seqfn, 0),
                        (wx.StaticText(self, -1, ""), 0, wx.EXPAND),
                        (self.cb_multilocus, 1, wx.EXPAND),
                        (wx.StaticText(self, -1, "Data Type"), 0, wx.ALIGN_RIGHT),
                        (self.datatype, 0),
                        (wx.StaticText(self, -1, "Initial Alignment"), 0, wx.ALIGN_RIGHT),
                        (self.checkbox_aligned, 0),
                        (self.tree_btn, 0, wx.LEFT|wx.EXPAND),
                        (self.txt_treefn, 0),
        ])

        self.ctrls.extend([self.seq_btn,
                           self.txt_seqfn,
                           self.tree_btn,
                           self.txt_treefn,
                           self.datatype])
        staticboxsizer.Add(sizer, 0, wx.CENTER, 0)
        self.Bind(wx.EVT_BUTTON, self.OnChooseSeq, self.seq_btn)
        self.Bind(wx.EVT_BUTTON, self.OnChooseTree, self.tree_btn)
        self.Bind(wx.EVT_COMBOBOX, self.OnDataType, self.datatype)
        self.Bind(wx.EVT_CHECKBOX, self.OnMultiLocus, self.cb_multilocus)
        return staticboxsizer

    def _create_workflow_settings_sizer(self):
        """
        returns a wx.StaticBoxSizer with the widgets that control pre and post
        processing of PASTA output.
        """
        staticboxsizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Workflow Settings"), wx.VERTICAL)
        sizer = wx.GridBagSizer(GRID_VGAP, GRID_HGAP)
        self.two_phase = wx.CheckBox(self, -1, "Two-Phase (not PASTA)")
        self.two_phase.Value = False
        self.raxml_after = wx.CheckBox(self, -1, "Extra RAxML Search")
        self.raxml_after.Value = False
        #self.trusted_data = wx.CheckBox(self, -1, "Trusted Data")
        #self.trusted_data.Value = True

        self.ctrls.extend([self.two_phase,
                           ])

        cr = 0
        sizer.Add(wx.StaticText(self, -1, "Algorithm"), (cr,0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.two_phase, (cr,1), flag=wx.EXPAND)

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Post-Processing"), (cr,0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.raxml_after, (cr,1), flag=wx.EXPAND)

        #cr += 1
        #sizer.Add(wx.StaticText(self, -1, "Input Validation"), (cr,0), flag=wx.ALIGN_LEFT )
        #sizer.Add(self.trusted_data, (cr,1), flag=wx.EXPAND)

        self.Bind(wx.EVT_CHECKBOX, self.OnTwoPhase, self.two_phase)

        staticboxsizer.Add(sizer, 0, wx.ALL, 0)
        return staticboxsizer
        
    def _create_pasta_settings_sizer(self):
        staticboxsizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "PASTA Settings"), wx.VERTICAL)
        sizer = wx.GridBagSizer(GRID_VGAP, GRID_HGAP)

#        preset_choices = ["SATe-II-fast", "SATe-II-ML", "SATe-II-simple", "(Custom)",]
#        self.cb_sate_presets = wx.ComboBox(self,
#                -1,
#                "SATe-II-ML",
#                choices=preset_choices,
#                style=wx.CB_READONLY)

        tree_and_alignment_choices = ["Final", "Best"]
        self.cb_tree_and_alignment = wx.ComboBox(self,
                -1,
                tree_and_alignment_choices[0],
                choices=tree_and_alignment_choices,
                style=wx.CB_READONLY)

        timelimit_list = list(map(str, [i/100.0 for i in range(1,10)] + [i/10.0 for i in range(1,10)] + list(range(1,73))))
        iterlimit_list = list(map(str, [1, 5, 10, 20, 50, 100, 200, 500, 1000]))
        self.rb_maxsub1 = wx.RadioButton(self, -1, "Percentage", name="frac", style=wx.RB_GROUP)
        self.rb_maxsub2 = wx.RadioButton(self, -1, "Size", name="size")
        self.cb_maxsub1 = wx.ComboBox(self, -1, "50", choices=list(map(str, list(range(1,51)))), style=wx.CB_READONLY)
        self.cb_maxsub2 = wx.ComboBox(self, -1, "200", choices=list(map(str, list(range(1,201)))), style=wx.CB_READONLY)

        self.ctrls.extend([self.rb_maxsub1,
                           self.cb_maxsub1,
                           self.rb_maxsub2,
                           self.cb_maxsub2
                           ])

        self.checkbox_stop_time = wx.CheckBox(self, -1, "Time Limit (hr)")
        self.checkbox_stop_iter = wx.CheckBox(self, -1, "Iteration Limit")
        self.cb_stop1 = wx.ComboBox(self, -1, "24", choices=timelimit_list, style=wx.CB_READONLY)
        self._iter_limits = [1, None]
        riv = RangedIntValidator(self._iter_limits[0], self._iter_limits[1])
        self.text_stop2 = wx.TextCtrl(self, -1, "8", validator=riv)
            
#        self.text_stop2.Bind(wx.EVT_KILL_FOCUS, lambda event : self.validate_iter_limit_text() and event.Skip())
#        self.blindmode = wx.CheckBox(self, -1, "Blind Mode Enabled")


#        self.blindmode.Value = True

        self.ctrls.extend([self.checkbox_stop_time,
                           self.cb_stop1,
                           self.checkbox_stop_iter,
                           self.text_stop2,
                           ])

        strategy_list = ["MinCluster","Centroid", "Longest"]
        self.cb_decomp = wx.ComboBox(self, -1, "MinCluster", choices=strategy_list, style=wx.CB_READONLY)

        self.ctrls.append(self.cb_decomp)
        self.pasta_settings_ctrl_list = []
        cr = 0

#        sizer.Add(wx.StaticText(self, -1, "Quick Set"), (cr, 0), flag=wx.ALIGN_LEFT )
#        sizer.Add(self.cb_sate_presets, (cr, 1), flag=wx.EXPAND)
#        self.pasta_settings_ctrl_list.append(self.cb_sate_presets)

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Max. Subproblem"), (cr,0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.rb_maxsub1, (cr,1), flag=wx.ALIGN_LEFT)
        sizer.Add(self.cb_maxsub1, (cr,2), flag=wx.EXPAND)
        self.pasta_settings_ctrl_list.extend([self.rb_maxsub1, self.cb_maxsub1])

        cr += 1
        sizer.Add(self.rb_maxsub2, (cr,1), flag=wx.ALIGN_LEFT)
        sizer.Add(self.cb_maxsub2, (cr,2), flag=wx.EXPAND)
        self.pasta_settings_ctrl_list.extend([self.rb_maxsub2, self.cb_maxsub2])

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Decomposition"), (cr,0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.cb_decomp, (cr,1), flag=wx.EXPAND)
        self.pasta_settings_ctrl_list.extend([self.cb_decomp])


#        cr += 1
#        sizer.Add(wx.StaticText(self, -1, "Apply Stop Rule"), (cr,0), flag=wx.ALIGN_LEFT )
#        sizer.Add(self.cb_apply_stop_rule, (cr,1), flag=wx.EXPAND)
#        self.pasta_settings_ctrl_list.extend([self.cb_apply_stop_rule])

#        cr += 1
#        sizer.Add(wx.StaticText(self, -1, "Stopping Rule"), (cr,0), flag=wx.ALIGN_LEFT )
#        sizer.Add(self.blindmode, (cr,1), flag=wx.EXPAND)
#        self.pasta_settings_ctrl_list.extend([self.blindmode])

        cr += 1
        sizer.Add(self.checkbox_stop_time, (cr,1), flag=wx.ALIGN_LEFT)
        sizer.Add(self.cb_stop1, (cr,2), flag=wx.EXPAND)
        self.pasta_settings_ctrl_list.extend([self.checkbox_stop_time, self.cb_stop1])

        cr += 1
        sizer.Add(self.checkbox_stop_iter, (cr,1), flag=wx.ALIGN_LEFT)
        sizer.Add(self.text_stop2, (cr,2), flag=wx.EXPAND)
        self.pasta_settings_ctrl_list.extend([self.checkbox_stop_iter, self.text_stop2])

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Return"), (cr, 0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.cb_tree_and_alignment, (cr, 1), flag=wx.EXPAND)
        self.pasta_settings_ctrl_list.extend([self.cb_tree_and_alignment])

        self.cb_maxsub1.Disable()
        self.cb_maxsub2.Disable()
        self.rb_maxsub1.Value = True
        self.cb_maxsub1.Enable()

        self.checkbox_stop_time.Value = False
        self.cb_stop1.Disable()
        self.text_stop2.Enable()
        self.checkbox_stop_iter.Value = True
        self.text_stop2.Value = "100"

        #self.Bind(wx.EVT_COMBOBOX, self.OnSatePresets, self.cb_sate_presets)
        #self.OnSatePresets(self.cb_sate_presets)

        #self.Bind(wx.EVT_CHECKBOX, self.OnBlindMode, self.blindmode)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnMaxSubproblem, self.rb_maxsub1)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnMaxSubproblem, self.rb_maxsub2)
        self.Bind(wx.EVT_CHECKBOX, self.OnTimeRuleCheckbox, self.checkbox_stop_time)
        self.Bind(wx.EVT_CHECKBOX, self.OnIterRuleCheckbox, self.checkbox_stop_iter)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_pasta_settings, self.cb_decomp)
        #self.Bind(wx.EVT_COMBOBOX, self._set_custom_pasta_settings, self.cb_apply_stop_rule)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_pasta_settings, self.cb_stop1)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_pasta_settings, self.text_stop2)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_pasta_settings, self.cb_tree_and_alignment)

        #cr += 1
        #presets = wx.ComboBox(self, -1, "1", choices=map(str, range(1,9)), style=wx.CB_READONLY)
        #sizer.Add(wx.StaticText(self, -1, "Preset Configuration"), (cr,0), flag=wx.ALIGN_LEFT )
        #sizer.Add(presets, (cr,1), flag=wx.EXPAND)

        staticboxsizer.Add(sizer, 0, wx.ALL, 0)
        return staticboxsizer

    def validate_iter_limit_text(self):
        field=self.text_stop2
        t = field.GetValue()
        if is_valid_int_str(t, self._iter_limits[0], self._iter_limits[1]):
            return True
        field.SetBackgroundColour("red")
        #wx.MessageBox(message='"Iteration Limit" must contain a positive integer',
        #                              caption='Input Error', style=wx.OK|wx.ICON_ERROR)
        field.SetFocus()
        field.Refresh()
        return False

    def _create_menu(self):
        self.menuBar = wx.MenuBar()
        self.menuFile = wx.Menu()
        self.menuHelp = wx.Menu()
        self.menuFileSaveLog = self.menuFile.Append(-1, "&Save Log...\tCtrl+S")
        self.menuFileExit = self.menuFile.Append(wx.ID_EXIT, "&Quit PASTA\tCtrl+Q")
        self.menuHelpHelp = self.menuHelp.Append( -1, "&Help")
        self.menuHelpAbout = self.menuHelp.Append(wx.ID_ABOUT, "&About PASTA")
        self.menuBar.Append(self.menuFile, "&File")
        self.menuBar.Append(self.menuHelp, "&Help")
        self.SetMenuBar(self.menuBar)
        self.Bind(wx.EVT_MENU, self.OnSaveLog, self.menuFileSaveLog)
        self.Bind(wx.EVT_MENU, self.OnExit, self.menuFileExit)
        self.Bind(wx.EVT_MENU, self.OnHelp, self.menuHelpHelp)
        self.Bind(wx.EVT_MENU, self.OnAbout, self.menuHelpAbout)
        self.Bind(wx.EVT_CLOSE, self.OnExit)

    def OnTreeEstimatorChange(self, event):
        self.set_char_model()

    def OnDataType(self, event):
        self.set_char_model()
    def set_char_model(self):
        if self.datatype.Value == "DNA" or self.datatype.Value == "RNA":
            self.cb_tools["model"].Clear()
            if self.cb_tools["treeestimator"].Value.lower() == "raxml":
                self.raxml_after.Value = False
                for model in self.raxml_dna_models:
                    self.cb_tools["model"].Append(model)
                    self.cb_tools["model"].SetStringSelection("GTRCAT")
            elif self.cb_tools["treeestimator"].Value.lower() == "fasttree":
                for model in self.fasttree_dna_models:
                    self.cb_tools["model"].Append(model)
                    self.cb_tools["model"].SetStringSelection("GTR+G20")
        elif self.datatype.Value == "Protein":
            self.cb_tools["model"].Clear()
            if self.cb_tools["treeestimator"].Value.lower() == "raxml":
                self.raxml_after.Value = False
                for model in self.raxml_prot_models:
                    self.cb_tools["model"].Append(model)
                    self.cb_tools["model"].SetStringSelection("PROTCATWAGF")
            elif self.cb_tools["treeestimator"].Value.lower() == "fasttree":
                for model in self.fasttree_prot_models:
                    self.cb_tools["model"].Append(model)
                    self.cb_tools["model"].SetStringSelection("JTT+G20")

    def OnSaveLog(self, event):
        dialog = wx.FileDialog(None, "Save Log", defaultFile=self.txt_jobname.Value, wildcard = "Log files (*.log)|*.log", style=wx.FD_OVERWRITE_PROMPT|wx.FD_SAVE)
        dialog.ShowModal()
        fn = dialog.GetPath()
        if len(fn) > 4:
            if not fn[-4:] == ".log":
                fn += ".log"
        else:
            fn += ".log"
        fc = open(fn, "w")
        fc.write(self.log.GetValue())
        fc.close()

    def OnMaxSubproblem(self, event):
        self._set_custom_pasta_settings(event)
        radio_selected = event.GetEventObject()
        if radio_selected.GetName() == "frac":
            self.cb_maxsub1.Enable()
            self.cb_maxsub2.Disable()
        elif radio_selected.GetName() == "size":
            self.cb_maxsub2.Enable()
            self.cb_maxsub1.Disable()
    def OnTwoPhase(self, event):
        """
        Called every time the 'Two-Phase' checkbox is clicked. The main action
            that needs to occur is the Disabling/Enabling of the PASTA settings 
            controls
        """
        if self.two_phase.Value:
            for c in self.pasta_settings_ctrl_list:
                c.Disable()
            self.cb_tools["merger"].Disable()
            self.tree_btn.Disable()
            self.txt_treefn.Disable()
            self.raxml_after.Disable()
        else:
            fragile_list = [self.cb_maxsub1, self.cb_maxsub2, self.cb_stop1, self.text_stop2]
            for c in self.pasta_settings_ctrl_list:
                if c not in fragile_list:
                    c.Enable()
            if self.rb_maxsub1.Value:
                self.cb_maxsub1.Enable()
            else:
                self.cb_maxsub2.Enable()
            if self.checkbox_stop_time.Value:
                self.cb_stop1.Enable()
            else:
                self.text_stop2.Enable()
            self.cb_tools["merger"].Enable()
            self.tree_btn.Enable()
            self.txt_treefn.Enable()
            self.raxml_after.Enable()


    def OnTimeRuleCheckbox(self, event):
        self._set_custom_pasta_settings(event)
        if self.checkbox_stop_time.Value:
            self.cb_stop1.Enable()
        else:
            self.cb_stop1.Disable()

    def OnIterRuleCheckbox(self, event):
        self._set_custom_pasta_settings(event)
        if self.checkbox_stop_iter.Value:
            self.text_stop2.Enable()
        else:
            self.text_stop2.Disable()

    def OnExit(self, event):
        if self.process is not None:
            wx.Process.Kill(self.pid, wx.SIGKILL)
        self._remove_config_file()
        self.Destroy()

    def OnHelp(self, event):
        import wx.html
        wx.FileSystem.AddHandler(wx.ZipFSHandler())
        def _addBook(filename):
            if not self.help.AddBook(filename, True):
                wx.MessageBox("Unable to open: " + filename, "Error", wx.OK|wx.ICON_EXCLAMATION)
        self.help = wx.html.HtmlHelpController(style = wx.html.HF_DEFAULT_STYLE^wx.html.HF_BOOKMARKS^wx.html.HF_INDEX)
        _addBook("help.zip")
        self.help.DisplayContents()

    def OnAbout(self, event):
        from wx.lib.wordwrap import wordwrap
        info = wx.AboutDialogInfo()
        info.SetName(PROGRAM_NAME)
        info.SetVersion(PROGRAM_VERSION)
        info.SetCopyright("Copyright (C) %s" % PROGRAM_YEAR)
        info.SetWebSite((PROGRAM_WEBSITE, "%s Homepage" % PROGRAM_NAME))
        info.SetLicense(PROGRAM_LICENSE)
        info.SetDescription(PROGRAM_DESCRIPTION)
        [info.AddDeveloper(i) for i in PROGRAM_AUTHOR]
        wx.AboutBox(info)

    def _show_error_dialog(self, error_msg, caption):
        """
        Puts up a modal dialog with a `error_msg` and `caption`
        destroys the dialog after the user clicks `OK`
        """
        error_msg_dlg = wx.MessageDialog(parent=self,  
                                     message=error_msg,
                                     caption=caption,
                                     style=wx.OK|wx.ICON_ERROR)
        error_msg_dlg.ShowModal()
        error_msg_dlg.Destroy()

    def OnChooseSeq(self, event):
        filepath = None
        parse_as_multilocus = self.cb_multilocus.Value
        if not parse_as_multilocus:
            dialog = wx.FileDialog(None, "Choose sequences...", wildcard = "FASTA files (*.fasta)|*.fasta|FASTA files (*.fas)|*.fas|FASTA files (*)|*", style=wx.FD_OPEN)
            dialog.ShowModal()
            self.txt_seqfn.SetValue( dialog.GetPath() )
            filepath = self._encode_arg(self.txt_seqfn.GetValue())
            if filepath and not self.txt_outputdir.GetValue():
                self.txt_outputdir.SetValue(os.path.dirname(os.path.abspath(filepath)))
        else:
            dialog = wx.DirDialog(None, "Choose directory for multiple sequence files", style=wx.FD_OPEN)
            dialog.ShowModal()
            self.txt_seqfn.SetValue( dialog.GetPath() )
            filepath = self._encode_arg(self.txt_seqfn.GetValue())
        if PARSING_FILES_IN_GUI and filepath:
            confirm_parse_dlg = wx.MessageDialog(parent=self,  
                                                 message="Do you want PASTA to read the data now? (this causes PASTA to customize some of the settings for your data).",
                                                 caption="Read input data now?", 
                                                 style=wx.OK|wx.CANCEL|wx.ICON_QUESTION)
            result = confirm_parse_dlg.ShowModal()
            confirm_parse_dlg.Destroy()
            if result == wx.ID_OK:
                progress_dialog = wx.ProgressDialog(title="Reading input data",
                                         message="Parsing data files        ",
                                         maximum=100,
                                         parent=self,
                                         style=wx.PD_AUTO_HIDE|wx.PD_APP_MODAL)
                progress_dialog.Update(1, "Beginning Parse")
                error_msg = None
                try:
                    if parse_as_multilocus:
                        fn_list = get_list_of_seq_filepaths_from_dir(filepath)
                    else:
                        fn_list = [filepath]
#                     if self.datatype.Value == "Protein":
#                         datatype_list = ["PROTEIN"]
#                     else:
                    datatype_list = ["DNA", "RNA", "PROTEIN"]
                    careful_parse = False
                    summary_stats = summary_stats_from_parse(fn_list, 
                                                             datatype_list,
                                                             None,
                                                             careful_parse=careful_parse)
                    progress_dialog.Update(100, "Done")
                except Exception as x:
                    try:
                        error_msg = "Problem reading the data:\n" + str(x.message)
                    except:
                        error_msg = "Unknown error encountered while reading the data."
                except:
                    error_msg = "Unknown error encountered while reading the data."
                if error_msg:
                    self._show_error_dialog(error_msg, caption="Input parsing error")
                    filepath = None
                    self._could_be_aligned = False
                    self.refresh_aligned_checkbox()
                else:
                    read_type = summary_stats[0]                    
                    if read_type == "PROTEIN":
                        self.datatype.SetValue("Protein")
                    else:
                        self.datatype.SetValue(read_type)
                    # Set defaults from "auto_defaults"
                    auto_opts = get_auto_defaults_from_summary_stats(summary_stats[0], summary_stats[1], summary_stats[2])

                    self._could_be_aligned = summary_stats[3]
                    self.refresh_aligned_checkbox()

                    auto_pasta_opts = auto_opts["sate"]
                    te_str = auto_pasta_opts["tree_estimator"].upper()
                    self.cb_tools["treeestimator"].SetStringSelection(te_str)
                    self.set_char_model()
                    if te_str == "FASTTREE":
                        te_opts = auto_opts['fasttree']
                    self.cb_tools["model"].SetStringSelection(te_opts["GUI_model"])
                    self.cb_tools["merger"].SetStringSelection(auto_pasta_opts["merger"].upper())
                    self.cb_tools["aligner"].SetStringSelection(auto_pasta_opts["aligner"].upper())
                    self.cb_ncpu.SetStringSelection(str(min(MAX_NUM_CPU, auto_pasta_opts["num_cpus"])))
                    
                    # Set max decomposition based on data set size (always move to actual # here)
                    self.rb_maxsub1.Value = False
                    self.cb_maxsub1.Disable()
                    self.rb_maxsub2.Value = True
                    self.cb_maxsub2.SetStringSelection(str(max(1, auto_pasta_opts["max_subproblem_size"])))
                    self.cb_maxsub2.Enable()
                    
                    bs = auto_pasta_opts["break_strategy"]
                    bs = {"mincluster":"MinCluster","centroid":"Centroid", "longest":"Longest"}[bs]
                    #bs = bs[0].upper() + bs[1:].lower()
                    self.cb_decomp.SetValue(bs)


                    self.cb_stop1.Disable()
                    self.checkbox_stop_iter.Value = True
                    self.cb_maxsub2.SetStringSelection(str(max(1, auto_pasta_opts["max_subproblem_size"])))
                    self.cb_maxsub2.Enable()
                    
                    if auto_pasta_opts['move_to_blind_on_worse_score']:
                        #self.blindmode.Value = True
                        t_l = auto_pasta_opts['after_blind_time_without_imp_limit']                        
                    else:
                        #self.blindmode.Value = False
                        t_l = auto_pasta_opts['time_limit']                        
                    
                    if t_l <= 0:
                        self.checkbox_stop_time.Value = False
                    else:
                        self.checkbox_stop_time.Value = True
                        
                        
#                     self.cb_apply_stop_rule.SetValue("After Last Improvement")
                    after_blind_it_lim = auto_pasta_opts['iter_limit']
                    self.text_stop2.SetValue(str(after_blind_it_lim))
                    
                    if self._could_be_aligned:
                        a_tag = "aligned"
                    else:
                        a_tag = "unaligned"

                    self.log.AppendText("Read %d file(s) with %s %s data. Total of %d taxa found.\n" % (len(fn_list), a_tag, read_type, summary_stats[2]))
                    by_file = summary_stats[1]
                    for n, fn in enumerate(fn_list):
                        t_c_tuple = by_file[n]
                        if self._could_be_aligned:
                            self.log.AppendText('  Parsing of the file "%s" returned %d sequences of length = %d\n' % (fn, t_c_tuple[0], t_c_tuple[1]))
                        else:
                            self.log.AppendText('  Parsing of the file "%s" returned %d sequences with longest length = %d\n' % (fn, t_c_tuple[0], t_c_tuple[1]))
                    
                progress_dialog.Destroy()
            else:
                self._could_be_aligned = True
                self.refresh_aligned_checkbox()
        if filepath:
            if not parse_as_multilocus:
                if filepath and not self.txt_outputdir.GetValue():
                    self.txt_outputdir.SetValue(os.path.dirname(os.path.abspath(filepath)))
            else:
                if filepath and not self.txt_outputdir.GetValue():
                    self.txt_outputdir.SetValue(os.path.abspath(filepath))
        else:
            self.txt_seqfn.SetValue("")
            
    def refresh_aligned_checkbox(self):
        self.checkbox_aligned.SetValue(self._could_be_aligned)
        if self._could_be_aligned:
            treefilename = self.txt_treefn.GetValue()
            if treefilename and os.path.isfile(treefilename):
                self.checkbox_aligned.Disable()
            else:
                self.checkbox_aligned.Enable()
        else:
            self.checkbox_aligned.Disable()


    def OnChooseTree(self, event):
        dialog = wx.FileDialog(None, "Choose tree...", wildcard = "Tree files (*.tre)|*.tre|Tree files (*.tree)|*.tree|Tree files (*.phy)|*.phy", style=wx.FD_OPEN)
        dialog.ShowModal()
        self.txt_treefn.SetValue( dialog.GetPath() )
        self.refresh_aligned_checkbox()

    def OnIdle(self, evt):
        if self.process is not None:
            stream = self.process.GetInputStream()
            if stream is not None and stream.CanRead():
                text = stream.read()
                self.log.AppendText(text)

            stream = self.process.GetErrorStream()
            if stream is not None and stream.CanRead():
                text = stream.read()
                self.log.AppendText(text)

    def OnProcessEnded(self, evt):
        stream = self.process.GetInputStream()
        if stream.CanRead():
            text = stream.read()
            self.log.AppendText(text)

        stream = self.process.GetErrorStream()
        if stream.CanRead():
            text = stream.read()
            self.log.AppendText(text)

        self.process.Destroy()
        self.process = None
        self.log.AppendText("Job %s is finished.\n" % self.txt_jobname.GetValue())
        self._remove_config_file()
        self._ReactivateOptions()
        self.statusbar.SetStatusText("PASTA Ready!")
        self.button.SetLabel("Start")

    def OnButton(self, event):
        if self.button.GetLabel() == "Start":
            self._OnStart()
        elif self.button.GetLabel() == "Stop":
            self._OnStop()
        else:
            raise ValueError("Button label %s not recognized.\n" % self.button.GetLabel() )

    def OnMultiLocus(self, event):
        if self.cb_multilocus.Value:
            self.seq_btn.SetLabel("Sequence files ...")
        else:
            self.seq_btn.SetLabel("Sequence file ...")
        self.txt_seqfn.SetValue("")

    def _FreezeOptions(self):
        self.prev_ctrls_status = []
        for ctrl in self.ctrls:
            self.prev_ctrls_status.append( ctrl.IsEnabled() )
            ctrl.Disable()

    def _ReactivateOptions(self):
        for i in range(len(self.ctrls)):
            self.ctrls[i].Enable(self.prev_ctrls_status[i])

    def _OnStart(self):


        if self.process is None:
            if (not self.checkbox_stop_time.Value) and (not self.checkbox_stop_iter.Value):
                self._show_error_dialog("Termination conditions are not set correctly. Either a time limit, an iteration limit, or both must be used.\n", caption="PASTA Settings Error")
                return
    
            if self.checkbox_stop_iter.Value and (not self.validate_iter_limit_text()):
                self._show_error_dialog("Iteration limit is not set correctly. Enter a positive integer in that field.\n", caption="PASTA Settings Error")
                return
            input_filename = self._encode_arg(self.txt_seqfn.GetValue())
            if not input_filename:
                self._show_error_dialog("Input sequence file(s) are required.\n", caption="PASTA Settings Error")
                return
            if not os.path.exists(input_filename):
                self._show_error_dialog('Input sequence file(s) are "%s" does not exist!.\n', caption="PASTA Settings Error")
                return
            if self.cb_multilocus.Value:
                if not os.path.isdir(input_filename):
                    self._show_error_dialog('Input sequence file specification should be a directory when multilocus model is used.\n', caption="PASTA Settings Error")
                    return
            elif not os.path.isfile(input_filename):
                self._show_error_dialog('Input sequence file must be a file when single-locus mode used.\n', caption="PASTA Settings Error")
                return

            cfg_success = self._create_config_file()
            if not cfg_success:
                return
            #command = [filemgr.quoted_file_path(x) for x in get_invoke_run_pasta_command()]
            command = get_invoke_run_pasta_command()
            treefilename = self._encode_arg(self.txt_treefn.GetValue())
            jobname = self._encode_arg(self.txt_jobname.GetValue())
            if not jobname:
                wx.MessageBox("Job name cannot be empty, it is REQUIRED by PASTA!", "WARNING", wx.OK|wx.ICON_WARNING)
                self._remove_config_file()
                return
            command.extend(["-i", filemgr.quoted_file_path(input_filename)])
            if treefilename and os.path.isfile(treefilename):
                command.extend(["-t", filemgr.quoted_file_path(treefilename)])
            command.extend(["-j", filemgr.quoted_file_path(jobname) ])
            if self.datatype.Value == "DNA":
                dt = "dna"
            elif self.datatype.Value == "RNA":
                dt = "rna"
            else:
                dt = "protein"
            command.extend(["-d", dt])
            command.extend(["%s" % filemgr.quoted_file_path(self.process_cfg_file)])
            if PASTA_GUI_ONLY_PRINTS_CONFIG:
                self.log.AppendText("Command is:\n  '%s'\n" % "' '".join(command))
                self.log.AppendText("config_file:\n#############################################################\n")
                for line in open(self.process_cfg_file, 'rU'):
                    self.log.AppendText(line)
                self.log.AppendText("#############################################################\n")
                self._remove_config_file()
                self.statusbar.SetStatusText("\n\nRun emulated!\n\n")
            else:
                print("running %s" %" ".join(command))
                self.process = wx.Process(self)
                self.process.Redirect()
                self.pid = wx.Execute( " ".join(command), wx.EXEC_ASYNC, self.process)
                self.button.SetLabel("Stop")
                self.statusbar.SetStatusText("PASTA Running!")
                self._FreezeOptions()

        else:
            self.log.AppendText("Job %s is still running!\n" % self.txt_jobname.GetValue())

    def _OnStop(self):
        if self.process is not None:
            self.log.AppendText("Job %s is terminated early.\n" % self.txt_jobname.GetValue())
            self.process.Kill(self.pid, wx.SIGKILL)
            self._remove_config_file()
            self._ReactivateOptions()
            self.button.SetLabel("Start")
            self.statusbar.SetStatusText("PASTA Ready!")
        else:
            self.log.AppendText("No active PASTA jobs to terminate!\n")

    def _encode_arg(self, arg, encoding='utf-8'):
        #if isinstance(arg, str):
        #    return arg.encode(encoding)
        return arg

    def _create_config_file(self):
        from pasta.configure import get_configuration
        cfg = get_configuration()

        #if self.txt_resultdir.Value:
        #    basefilename = os.path.basename(self.txt_seqfn.GetValue())
        #    jobname = self.txt_jobname.GetValue()
        #    resultdir = self.txt_resultdir.Value
        #    cfg.commandline.output = os.path.join(resultdir, basefilename+"_%s.aln" % jobname )
        #    cfg.commandline.result = os.path.join(resultdir, basefilename+"_%s.tre" % jobname )

        cfg.sate.aligner = self.cb_tools["aligner"].Value
        cfg.sate.tree_estimator = self.cb_tools["treeestimator"].Value
        if self.cb_tools["treeestimator"].Value.lower() == "raxml":
            cfg.raxml.model = self.cb_tools["model"].Value
        else:
            model_desc = self.cb_tools["model"].Value
            if model_desc == "GTR+G20":
                cfg.fasttree.model = "-gtr -gamma"
            elif model_desc == "GTR+CAT":
                cfg.fasttree.model = "-gtr"
            elif model_desc == "JC+G20":
                cfg.fasttree.model = "-gamma"
            elif model_desc == "JC+CAT":
                cfg.fasttree.model = ""
            elif model_desc == "JTT+G20":
                cfg.fasttree.model = "-gamma"
            elif model_desc == "JTT+CAT":
                cfg.fasttree.model = ""
            elif model_desc == "WAG+G20":
                cfg.fasttree.model = "-wag -gamma"
            elif model_desc == "WAG+CAT":
                cfg.fasttree.model = "-wag"
            else:
                raise Exception("Unrecognized model: %s" % model_desc)
        cfg.commandline.keeptemp = True
        cfg.commandline.keepalignmenttemps = True
        if self.checkbox_aligned.Value:
            cfg.commandline.aligned = True
        #cfg.commandline.untrusted = not bool(self.trusted_data.Value)
        if self.cb_multilocus.Value:
            cfg.commandline.multilocus = True
        
        if self.two_phase.Value:
            cfg.commandline.two_phase = True
            cfg.commandline.raxml_search_after = False
        else:
            cfg.commandline.two_phase = False
            cfg.commandline.raxml_search_after = bool(self.raxml_after.Value)
            cfg.sate.merger = self.cb_tools["merger"].Value
            cfg.sate.break_strategy = self.cb_decomp.Value
            cfg.sate.start_tree_search_from_current = True
            if self.rb_maxsub1.Value:
                cfg.sate.max_subproblem_frac = float(self.cb_maxsub1.Value)/100.0
                cfg.sate.max_subproblem_size = 0
            elif self.rb_maxsub2.Value:
                cfg.sate.max_subproblem_size = self.cb_maxsub2.Value
                cfg.sate.max_subproblem_frac = 0.0
    
    
            cfg.sate.time_limit = -1
            cfg.sate.iter_limit = -1
            cfg.sate.after_blind_time_without_imp_limit = -1
            cfg.sate.after_blind_iter_without_imp_limit = -1
#             if True:
#                 cfg.pasta.move_to_blind_on_worse_score = True
#                 if self.cb_apply_stop_rule.GetValue() == "After Last Improvement":
#                     if self.checkbox_stop_time.Value:
#                         cfg.pasta.after_blind_time_without_imp_limit = float(self.cb_stop1.Value)*3600
#                     if self.checkbox_stop_iter.Value:
#                         cfg.pasta.after_blind_iter_without_imp_limit = int(self.text_stop2.Value)
#                 else:
#                     if self.checkbox_stop_time.Value:
#                         cfg.pasta.time_limit = float(self.cb_stop1.Value)*3600
#                     if self.checkbox_stop_iter.Value:
#                         cfg.pasta.iter_limit = int(self.text_stop2.Value)
#             else:
            if self.checkbox_stop_time.Value:
                cfg.sate.time_limit = float(self.cb_stop1.Value)*3600
            if self.checkbox_stop_iter.Value:
                cfg.sate.iter_limit = int(self.text_stop2.Value)
            cfg.sate.return_final_tree_and_alignment = self.cb_tree_and_alignment.GetValue() == "Final"

        cfg.sate.output_directory = self._encode_arg(self.txt_outputdir.GetValue())
        cfg.sate.num_cpus = self.cb_ncpu.Value
        max_mb = self.txt_maxmb.GetValue()
        if not self.validate_max_mb(max_mb):
            return False
        cfg.sate.max_mem_mb = max_mb

        # this creates a file that cannot be deleted while the Python
        # process is running (under the mess that is called "Windows")
        #tf, self.process_cfg_file = tempfile.mkstemp(dir=pasta_home_dir(),
        #        suffix="_internal.cfg")

        tf = tempfile.NamedTemporaryFile(suffix="_internal.cfg", dir=pasta_home_dir())
        self.process_cfg_file = tf.name
        tf.close()
        cfg.save_to_filepath(self.process_cfg_file)
        return True

    def _remove_config_file(self):
        if "PASTA_GUIDEVMODE" in os.environ:
            return
        if self.process_cfg_file and os.path.exists(self.process_cfg_file):
            try:
                os.remove(self.process_cfg_file)
            except:
                # on windows, the config file sometimes cannot be deleted
                # ("...because it is being used by another process...")
                # resulting in an exception being thrown.
                # so we just:
                pass

class PastaApp(wx.PySimpleApp):
    def OnInit(self):
        self.frame = PastaFrame(size=wx.Display().GetClientArea())
        self.frame.Show(True)
        self.SetTopWindow(self.frame)
        return True

def main_gui():
    app = PastaApp()
    app.MainLoop()

ICO_STR = """AAABAAMAEBAAAAAAIABoBAAANgAAACAgAAAAACAAqBAAAJ6EAAAwMAAAAAAgAKglAABGFQAAKAAA\nABAAAAAgAAAAAQAgAAAAAABABAAAAAAAAAAAAAAAAAAAAAAAAAAAAGsAAADvAAAAqQAAAEUAAACX\n////Af///wEAAAC3AAAAJQAAAGsAAABv////Af///wEAAACHAAAA8QAAAGkAAAD1AAAA/wAAAP8A\nAABpAAAA4f///wEAAAALAAAA/wAAABkAAACPAAAAk////wEAAAAVAAAA+wAAAP8AAADXAAAA8wAA\nALsAAAD7AAAAgQAAAPsAAAAlAAAAQwAAAPkAAAADAAAAjwAAAJP///8BAAAAUQAAAO0AAABfAAAA\ntwAAAGv///8BAAAAsQAAAHsAAAD/AAAA/wAAAP8AAADd////AQAAAI8AAACT////AQAAAHUAAACl\n////AQAAAB3///8B////AQAAAKcAAAB7AAAA6QAAAP8AAAD/AAAAv////wEAAACPAAAAk////wEA\nAACHAAAAsQAAAFMAAABT////AQAAABcAAADlAAAAcwAAAM0AAACvAAAAwQAAAKP///8BAAAAjwAA\nAJP///8BAAAAjwAAAP8AAAD/AAAA/wAAAB0AAADfAAAA/wAAAFsAAACvAAAAawAAAJMAAACH////\nAQAAAI8AAACT////AQAAAIsAAADnAAAAzwAAAPkAAACZAAAA/wAAAP0AAAAhAAAAkwAAAIUAAACv\nAAAAaf///wEAAACPAAAAk////wEAAAB7AAAAmQAAACMAAADrAAAA2QAAAP8AAACF////AQAAAHUA\nAAChAAAAyQAAAE3///8BAAAAjwAAAJP///8BAAAAWwAAANUAAABjAAAAzQAAAPMAAABv////Af//\n/wEAAABZAAAAvQAAAOUAAAAv////AQAAAI8AAACT////AQAAACMAAAD/AAAA/wAAAJUAAAD9AAAA\nE////wH///8BAAAAOwAAANsAAAD7AAAAEf///wEAAACPAAAAk////wH///8BAAAAtwAAAP0AAAAz\nAAAA9QAAACsAAAA3AAAAKQAAAB8AAAD9AAAA8wAAAAMAAAA3AAAApwAAAKkAAAA3AAAAAwAAAAsA\nAAAp////AQAAANkAAADvAAAA/QAAADMAAAAFAAAA+wAAANcAAAAJAAAA/wAAAP8AAAD/AAAA/wAA\nAAsAAAAFAAAAXf///wEAAACdAAAA/wAAAP8AAAAz////AQAAAOMAAAC5AAAACQAAAP8AAAD/AAAA\n/wAAAP8AAAAL////AQAAAKf///8BAAAAJQAAAMUAAACPAAAAC////wEAAAB3AAAAXwAAAAUAAACV\nAAAAlQAAAJUAAACVAAAAB////wEAAACZAAAAIf///wH///8B////Af///wH///8B////Af///wH/\n//8B////Af///wH///8B////Af///wH///8BAAAAXQAAAGsAAP//AAD//wAA//8AAP//AAD//wAA\n//8AAP//AAD//wAA//8AAP//AAD//wAA//8AAP//AAD//wAA//8AAP//KAAAACAAAABAAAAAAQAg\nAAAAAACAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAkAAABTAAAAyQAAAPkAAAC7AAAAMwAAAAcAAABb\nAAAAgQAAAEX///8B////Af///wH///8BAAAAbQAAAIEAAAA3////AQAAAB0AAABzAAAAdQAAAB//\n//8B////Af///wH///8BAAAAFwAAAJUAAAD3AAAA0QAAAFUAAAAJAAAAZwAAAOcAAAD/AAAA/wAA\nAP0AAAC1AAAABQAAAK0AAAD/AAAAmf///wH///8B////Af///wEAAADrAAAA/wAAAF3///8BAAAA\nOwAAAOUAAADnAAAAP////wH///8B////Af///wEAAAB5AAAA9wAAAP8AAAD/AAAA5wAAAF8AAADp\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAxAAAAkQAAAP8AAAC3////Af///wH///8BAAAACwAAAP8A\nAAD/AAAAPf///wEAAAA7AAAA5QAAAOcAAAA/////Af///wH///8BAAAADQAAAO0AAAD/AAAA/wAA\nAP8AAAD/AAAArwAAAOkAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAHEAAAB3AAAA/wAAAM////8B////\nAf///wEAAAAjAAAA/wAAAP8AAAAj////AQAAADsAAADlAAAA5wAAAD////8B////Af///wEAAABF\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAACvAAAA6QAAAP8AAAD/AAAA/wAAAP8AAAD/AAAApwAAAFkA\nAAD/AAAA7////wH///8B////AQAAAEEAAAD/AAAA+wAAAAf///8BAAAAOwAAAOUAAADnAAAAP///\n/wH///8B////AQAAAI0AAAD/AAAA+wAAAKsAAAC3AAAA9QAAAK8AAADpAAAA+QAAAIUAAABpAAAA\n8QAAAP8AAAC5AAAASwAAAP8AAAD9AAAASQAAAEcAAABHAAAAgwAAAP8AAADp////Af///wEAAAA7\nAAAA5QAAAOcAAAA/////Af///wEAAAAHAAAArQAAAP8AAAC/AAAACQAAABMAAACPAAAAqwAAANUA\nAABX////Af///wEAAAB7AAAA/wAAAMUAAAA3AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAMn///8B////AQAAADsAAADlAAAA5wAAAD////8B////AQAAABsAAADFAAAA/QAAAGP///8B////\nAQAAABMAAABXAAAAeQAAAAv///8B////AQAAAFMAAAD5AAAAywAAACcAAAD7AAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAArf///wH///8BAAAAOwAAAOUAAADnAAAAP////wH///8BAAAAJQAAAM8A\nAADvAAAARf///wH///8B////AQAAAAkAAAAF////Af///wH///8BAAAATQAAAPcAAADPAAAAJwAA\nAOEAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACN////Af///wEAAAA7AAAA5QAAAOcAAAA/////\nAf///wEAAAAxAAAA2wAAAOMAAAA5////Af///wH///8B////Af///wH///8B////Af///wEAAABd\nAAAA/QAAANEAAAAnAAAAxwAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAHP///8B////AQAAADsA\nAADlAAAA5wAAAD////8B////AQAAADUAAADfAAAA8wAAALkAAACnAAAApwAAAKcAAACl////Af//\n/wH///8BAAAAAwAAAJ8AAAD/AAAAywAAACEAAACnAAAA/wAAAPEAAADbAAAA3wAAAPkAAAD9AAAA\nVf///wH///8BAAAAOwAAAOUAAADnAAAAP////wH///8BAAAAOQAAAOMAAAD/AAAA/wAAAP8AAAD/\nAAAA/wAAAP////8B////AQAAAAMAAABTAAAA9QAAAP8AAADFAAAAGwAAAI0AAAD/AAAAuwAAADMA\nAABTAAAA3QAAAPEAAABJ////Af///wEAAAA7AAAA5QAAAOcAAAA/////Af///wEAAAA5AAAA4wAA\nAP8AAAD/AAAA/wAAAP8AAAD/AAAA/////wEAAAAHAAAAjwAAAPUAAAD/AAAA/wAAALkAAAARAAAA\nawAAAP8AAAC5AAAADwAAADkAAADhAAAA4QAAADf///8B////AQAAADsAAADlAAAA5wAAAD////8B\n////AQAAADcAAADhAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD3AAAAAwAAAGkAAAD5AAAA/wAAAP8A\nAAD/AAAAoQAAAAMAAABVAAAA+wAAAMUAAAAbAAAARQAAAO0AAADVAAAAK////wH///8BAAAAOwAA\nAOUAAADnAAAAP////wH///8BAAAAMwAAAN0AAADxAAAAqwAAAJUAAACrAAAA/wAAAPEAAAAvAAAA\n1wAAAP8AAAD/AAAA/wAAAP8AAABn////AQAAAEUAAADtAAAA1QAAACsAAABVAAAA+wAAAMMAAAAb\n////Af///wEAAAA7AAAA5QAAAOcAAAA/////Af///wEAAAArAAAA0wAAAOcAAAA/////AQAAADsA\nAAD/AAAA3wAAAGUAAAD7AAAA/wAAAP8AAAD/AAAA9wAAAB3///8BAAAANwAAAOEAAADhAAAANwAA\nAGkAAAD/AAAAtwAAAA3///8B////AQAAADsAAADlAAAA5wAAAD////8B////AQAAACEAAADLAAAA\n9QAAAEv///8BAAAATwAAAP8AAADPAAAAnwAAAP8AAAD/AAAA/wAAAPUAAAB/////Af///wEAAAAn\nAAAA0QAAAPEAAABHAAAAhwAAAP8AAACj////Af///wH///8BAAAAOwAAAOUAAADnAAAAP////wH/\n//8BAAAADwAAALkAAAD/AAAAfwAAAAMAAACFAAAA/wAAAKsAAADLAAAA/wAAAP8AAAD/AAAAiQAA\nABP///8B////AQAAABsAAADDAAAA+wAAAFMAAACfAAAA/wAAAIv///8B////Af///wEAAAA7AAAA\n5QAAAOcAAAA/////Af///wEAAAADAAAAowAAAP8AAADXAAAAPQAAAMcAAAD/AAAAhwAAAOMAAAD/\nAAAA4wAAAFv///8B////Af///wH///8BAAAACQAAALMAAAD/AAAAbQAAAL0AAAD/AAAAa////wH/\n//8B////AQAAADsAAADlAAAA5wAAAD////8B////Af///wEAAABlAAAA/wAAAP8AAAD/AAAA/wAA\nAPcAAABTAAAA7QAAAP8AAAB3AAAAB////wH///8B////Af///wH///8BAAAAowAAAP8AAACHAAAA\n1QAAAP8AAABR////Af///wH///8BAAAAOwAAAOUAAADnAAAAP////wH///8B////AQAAACUAAAD/\nAAAA/wAAAP8AAAD/AAAA2wAAADEAAAD7AAAA/wAAAC////8B////Af///wH///8B////Af///wEA\nAACDAAAA/wAAAKUAAADzAAAA/wAAAC////8B////Af///wEAAAA7AAAA5QAAAOcAAAA/////Af//\n/wH///8B////AQAAALEAAAD/AAAA/wAAAP8AAACZAAAABwAAAPsAAAD/AAAAGf///wH///8B////\nAf///wH///8B////AQAAAGkAAAD/AAAAywAAAP8AAAD/AAAAFf///wH///8B////AQAAADsAAADl\nAAAA5wAAAD////8B////Af///wH///8BAAAARwAAAOUAAAD/AAAA8wAAAC3///8BAAAA9QAAAP8A\nAAAx////Af///wEAAAAdAAAAPf///wH///8BAAAASQAAAP8AAAD7AAAA/wAAAPX///8B////Af//\n/wH///8BAAAAOwAAAOUAAADnAAAAP////wH///8B////Af///wH///8BAAAAKQAAAHMAAAAx////\nAf///wEAAADlAAAA/wAAAHUAAAADAAAAGwAAAKUAAABn////Af///wEAAAAvAAAA/wAAAP8AAAD/\nAAAA2////wEAAAAHAAAAawAAAGsAAACNAAAA7wAAAPEAAACPAAAAawAAAGsAAAAJ////Af///wH/\n//8B////Af///wH///8B////AQAAAMcAAAD/AAAA9QAAAMkAAAD1AAAA/wAAAGf///8B////AQAA\nAA8AAAD/AAAA/wAAAP8AAAC7////AQAAAA8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAABf///8B////AQAAAAkAAABDAAAAD////wH///8BAAAAowAAAP8AAAD/AAAA/wAAAP8AAAD/\nAAAAZ////wH///8B////AQAAAPMAAAD/AAAA/wAAAKH///8BAAAADwAAAP8AAAD/AAAA/wAAAP8A\nAAD/AAAA/wAAAP8AAAD/AAAAF////wH///8BAAAACwAAAN8AAABB////Af///wEAAABlAAAA+wAA\nAP8AAAD/AAAA/wAAAP8AAABn////Af///wH///8BAAAA0wAAAP8AAAD/AAAAf////wEAAAAPAAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAX////Af///wH///8BAAAAwwAAAIP///8B\n////AQAAADcAAADfAAAA/wAAAP8AAAD/AAAA/wAAAGf///8B////Af///wEAAAC5AAAA/wAAAP8A\nAABl////AQAAAA8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAABf///8B////Af//\n/wEAAACPAAAAyf///wH///8BAAAABQAAAH0AAAD9AAAA/wAAAPkAAADFAAAAKf///wH///8B////\nAQAAAI0AAADrAAAA5QAAAEn///8BAAAADwAAAOsAAADrAAAA6wAAAOsAAADrAAAA6wAAAOsAAADr\nAAAAFf///wH///8B////AQAAAFUAAAD5AAAAIf///wH///8BAAAAEQAAAHsAAACbAAAAYQAAAB//\n//8B////Af///wH///8BAAAAJQAAAEEAAAA9AAAAE////wEAAAAFAAAAQQAAAEEAAABBAAAAQQAA\nAEEAAABBAAAAQQAAAEEAAAAH////Af///wH///8BAAAANwAAAOEAAABj////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wEAAAAXAAAAwQAAAK8A\nAAAL////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAQAAAAMAAACXAAAAxwAAACcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACgAAAAwAAAAYAAAAAEAIAAAAAAAgCUAAAAA\nAAAAAAAAAAAAAAAAAAD///8B////AQAAAC0AAACTAAAA4QAAAP8AAADNAAAAOf///wH///8BAAAA\nCwAAAEEAAABBAAAAQQAAABH///8B////Af///wH///8B////Af///wEAAAAvAAAAQQAAAEEAAAAt\n////Af///wH///8BAAAAKwAAAEEAAABBAAAAL////wH///8B////Af///wH///8B////Af///wH/\n//8BAAAACwAAAJsAAAD3AAAA7QAAAJ8AAAAx////Af///wH///8BAAAARwAAAPUAAAD/AAAA/wAA\nAP8AAAD/AAAA9QAAADH///8BAAAAHQAAAP8AAAD/AAAA/wAAAFP///8B////Af///wH///8B////\nAf///wEAAADPAAAA/wAAAP8AAACf////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B\n////Af///wH///8B////Af///wH///8BAAAAtwAAAP8AAAD/AAAA/wAAAP8AAAD1AAAARf///wEA\nAABZAAAA+wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAMH///8BAAAAAwAAAPsAAAD/AAAA/wAA\nAG////8B////Af///wH///8B////Af///wEAAADrAAAA/wAAAP8AAACD////Af///wH///8BAAAA\nrwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////Af///wEAAABjAAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAA+QAAAEcAAADfAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8A\nAAAv////AQAAAOEAAAD/AAAA/wAAAIn///8B////Af///wH///8B////AQAAAAcAAAD/AAAA/wAA\nAP8AAABl////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////\nAf///wEAAADVAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAIcAAADfAAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAB7////AQAAAMUAAAD/AAAA/wAAAKX///8B////Af///wH/\n//8B////AQAAACEAAAD/AAAA/wAAAP8AAABJ////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf//\n/wH///8B////Af///wH///8B////AQAAADkAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAAIcAAADfAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAC/////AQAAAKcAAAD/\nAAAA/wAAAMH///8B////Af///wH///8B////AQAAAD0AAAD/AAAA/wAAAP8AAAAr////Af///wH/\n//8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////AQAAAH8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAIcAAADfAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAAP8AAADv////AQAAAIsAAAD/AAAA/wAAAN3///8B////Af///wH///8B////AQAAAFkAAAD/\nAAAA/wAAAP8AAAAN////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH/\n//8B////AQAAAMMAAAD/AAAA/wAAAP8AAAD9AAAA5QAAAP8AAAD/AAAA/wAAAIcAAADfAAAA/wAA\nAP8AAAD/AAAA/QAAAP8AAAD/AAAA/wAAAP8AAAD/AAAAHQAAAG0AAAD/AAAA/wAAAPf///8B////\nAf///wH///8B////AQAAAHUAAAD/AAAA/wAAAPH///8B////Af///wH///8BAAAArwAAAP8AAAD/\nAAAAuf///wH///8B////Af///wH///8B////AQAAAPMAAAD/AAAA/wAAAO8AAAAhAAAAAwAAAFMA\nAADVAAAA/wAAAIcAAADfAAAA/wAAAPMAAABlAAAACwAAACcAAADhAAAA/wAAAP8AAAD/AAAAOQAA\nAE8AAAD/AAAA/wAAAP8AAABxAAAAaQAAAGkAAABpAAAAaQAAALkAAAD/AAAA/wAAANP///8B////\nAf///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAHQAAAP8AAAD/\nAAAA/wAAAHX///8B////Af///wEAAAAfAAAA7wAAAIcAAADfAAAA9QAAADX///8B////Af///wEA\nAABTAAAA/wAAAP8AAAD/AAAASwAAADMAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAALf///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////\nAf///wH///8BAAAARwAAAP8AAAD/AAAA/wAAACX///8B////Af///wH///8BAAAAUQAAAIcAAADf\nAAAAVf///wH///8B////Af///wEAAAALAAAA+wAAAP8AAAD/AAAAXwAAABUAAAD/AAAA/wAAAP8A\nAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAJn///8B////Af///wH///8BAAAArwAA\nAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAYQAAAP8AAAD/AAAA5////wH///8B////\nAf///wH///8B////AQAAACkAAACL////Af///wH///8B////Af///wH///8BAAAA6QAAAP8AAAD/\nAAAAZ////wEAAAD3AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAHv/\n//8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAdwAA\nAP8AAAD/AAAAw////wH///8B////Af///wH///8B////Af///wEAAAAJ////Af///wH///8B////\nAf///wH///8BAAAA3QAAAP8AAAD/AAAAb////wEAAADbAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAF////8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH/\n//8B////Af///wH///8BAAAAiwAAAP8AAAD/AAAAr////wH///8B////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8BAAAA8wAAAP8AAAD/AAAAd////wEAAAC9AAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAEH///8B////Af///wH///8B\nAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAmwAAAP8AAAD/AAAAm////wH/\n//8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wEAAAAVAAAA/wAA\nAP8AAAD/AAAAc////wEAAAChAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAACX///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B\nAAAAoQAAAP8AAAD/AAAA/QAAAPkAAAD5AAAA+QAAAPkAAAD5AAAA+QAAAPf///8B////Af///wH/\n//8B////Af///wEAAAA/AAAA/wAAAP8AAAD/AAAAZf///wEAAACDAAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAAf///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAA\nuf///wH///8B////Af///wH///8BAAAApwAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/\nAAAA/wAAAP////8B////Af///wH///8B////AQAAAA0AAADTAAAA/wAAAP8AAAD/AAAAV////wEA\nAABnAAAA/wAAAP8AAAD7AAAAkQAAAJEAAACRAAAAwQAAAP8AAAD/AAAA6f///wH///8B////Af//\n/wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAArQAAAP8AAAD/AAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP////8B////Af///wH///8BAAAABwAAALEAAAD/\nAAAA/wAAAP8AAAD/AAAAS////wEAAABJAAAA/wAAAP8AAAD/AAAACf///wH///8BAAAAgwAAAP8A\nAAD/AAAAzf///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af//\n/wH///8BAAAAqwAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP3///8B////\nAf///wEAAAAlAAAA0wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAAO////wEAAAArAAAA/wAAAP8AAAD/\nAAAAI////wH///8BAAAAnwAAAP8AAAD/AAAAr////wH///8B////Af///wH///8BAAAArwAAAP8A\nAAD/AAAAuf///wH///8B////Af///wH///8BAAAApQAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAPX///8B////AQAAAB8AAADjAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\nFf///wEAAAAPAAAA/wAAAP8AAAD/AAAAP////wH///8BAAAAuQAAAP8AAAD/AAAAk////wH///8B\n////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAnwAAAP8A\nAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAO////8BAAAABwAAANsAAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAADj////Af///wH///8BAAAA8QAAAP8AAAD/AAAAW////wH///8BAAAA\n1QAAAP8AAAD/AAAAdf///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B\n////Af///wH///8BAAAAlwAAAP8AAAD/AAAAwwAAAGEAAABhAAAAYQAAAI8AAAD/AAAA/wAAAOf/\n//8BAAAAYwAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACx////Af///wH///8BAAAA1QAA\nAP8AAAD/AAAAd////wH///8BAAAA8QAAAP8AAAD/AAAAV////wH///8B////Af///wH///8BAAAA\nrwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAhQAAAP8AAAD/AAAAsf///wH///8B\n////AQAAAFUAAAD/AAAA/wAAANf///8BAAAA1wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8A\nAABv////Af///wH///8BAAAAtwAAAP8AAAD/AAAAkf///wEAAAALAAAA/wAAAP8AAAD/AAAAO///\n/wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAA\nbwAAAP8AAAD/AAAAyf///wH///8B////AQAAAGEAAAD/AAAA/wAAAMMAAAAtAAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAOsAAAAL////Af///wH///8BAAAAmQAAAP8AAAD/AAAArf///wEA\nAAAnAAAA/wAAAP8AAAD/AAAAHf///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf//\n/wH///8B////Af///wH///8BAAAAWQAAAP8AAAD/AAAA6////wH///8B////AQAAAH8AAAD/AAAA\n/wAAAK8AAABfAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAHf///8B////Af///wH///8B\nAAAAfQAAAP8AAAD/AAAAyf///wEAAABBAAAA/wAAAP8AAAD9AAAAA////wH///8B////Af///wH/\n//8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAOQAAAP8AAAD/AAAA/wAA\nACX///8B////AQAAALEAAAD/AAAA/wAAAI8AAACRAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\nnwAAAAP///8B////Af///wH///8BAAAAXwAAAP8AAAD/AAAA5f///wEAAABdAAAA/wAAAP8AAADj\n////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH/\n//8BAAAAEQAAAP8AAAD/AAAA/wAAAHf///8BAAAACwAAAO8AAAD/AAAA/wAAAGcAAAC/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAACN////Af///wH///8B////Af///wH///8BAAAAQwAAAP8AAAD/AAAA\n/QAAAAMAAAB3AAAA/wAAAP8AAADF////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/\nAAAAuf///wH///8B////Af///wH///8B////AQAAAOsAAAD/AAAA/wAAAOsAAAA/AAAAkQAAAP8A\nAAD/AAAA/wAAAD8AAADRAAAA/wAAAP8AAAD/AAAA7QAAAEn///8B////Af///wH///8B////Af//\n/wH///8BAAAAJQAAAP8AAAD/AAAA/wAAABsAAACTAAAA/wAAAP8AAACp////Af///wH///8B////\nAf///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////AQAAAK8AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA+wAAAAsAAADdAAAA/wAAAP8AAAD3AAAAO////wH/\n//8B////Af///wH///8B////Af///wH///8BAAAACQAAAP8AAAD/AAAA/wAAADcAAACvAAAA/wAA\nAP8AAACL////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////\nAf///wH///8B////AQAAAGsAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAAxf///wEAAADp\nAAAA/wAAAP8AAACB////Af///wH///8B////Af///wH///8B////Af///wH///8B////AQAAAOsA\nAAD/AAAA/wAAAFMAAADJAAAA/wAAAP8AAABv////Af///wH///8B////Af///wH///8BAAAArwAA\nAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////AQAAAB8AAAD9AAAA/wAAAP8AAAD/AAAA\n/wAAAP8AAAD/AAAAef///wEAAAD3AAAA/wAAAP8AAABT////Af///wH///8B////Af///wH///8B\n////Af///wH///8B////AQAAAM0AAAD/AAAA/wAAAG8AAADlAAAA/wAAAP8AAABR////Af///wH/\n//8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////Af//\n/wEAAACzAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD7AAAAGf///wEAAAD9AAAA/wAAAP8AAAAt////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////AQAAAK8AAAD/AAAA/wAAAI0AAAD9\nAAAA/wAAAP8AAAAz////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH/\n//8B////Af///wH///8B////Af///wEAAAA5AAAA/QAAAP8AAAD/AAAA/wAAAP8AAACf////Af//\n/wEAAAD3AAAA/wAAAP8AAAAj////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAQAAAJMAAAD/AAAA/wAAAMEAAAD/AAAA/wAAAP8AAAAX////Af///wH///8B////Af///wH///8B\nAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////Af///wH///8BAAAAiQAAAP8A\nAAD/AAAA/wAAAOUAAAAV////Af///wEAAADzAAAA/wAAAP8AAAA3////Af///wH///8B////Af//\n/wEAAAA7////Af///wH///8B////AQAAAHUAAAD/AAAA/wAAAPUAAAD/AAAA/wAAAPn///8B////\nAf///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B\n////Af///wH///8B////AQAAAF0AAACzAAAAnQAAAB3///8B////Af///wEAAADnAAAA/wAAAP8A\nAABr////Af///wH///8B////AQAAAHsAAACZ////Af///wH///8B////AQAAAFkAAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAN3///8B////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAA\nuf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wEAAADPAAAA/wAAAP8AAADPAAAACf///wEAAAAFAAAAbQAAAP8AAACZ////Af///wH/\n//8B////AQAAADsAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAL////8B////AQAAAA0AAAChAAAAoQAA\nAKEAAAChAAAA4QAAAP8AAAD/AAAA5QAAAKEAAAChAAAAoQAAAKEAAAAV////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wEAAAC3AAAA/wAAAP8AAAD/AAAA0wAAAJ0AAADp\nAAAA/wAAAP8AAACZ////Af///wH///8B////AQAAAB0AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAKH/\n//8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAAh////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wEAAACTAAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACZ////Af///wH///8B////AQAAAAUAAAD9\nAAAA/wAAAP8AAAD/AAAA/wAAAIX///8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8A\nAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAh////Af///wH///8B////AQAAACkAAADpAAAAhf//\n/wH///8B////Af///wEAAABlAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACZ////\nAf///wH///8B////Af///wEAAADjAAAA/wAAAP8AAAD/AAAA/wAAAGf///8B////AQAAABcAAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAh////Af///wH/\n//8B////AQAAAAMAAADxAAAA3f///wH///8B////Af///wEAAAAtAAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAACZ////Af///wH///8B////Af///wEAAADHAAAA/wAAAP8AAAD/AAAA\n/wAAAEn///8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAAh////Af///wH///8B////Af///wEAAAC5AAAA/wAAACv///8B////Af///wH/\n//8BAAAA4QAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACZ////Af///wH///8B////Af//\n/wEAAACpAAAA/wAAAP8AAAD/AAAA/wAAAC3///8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAh////Af///wH///8B////Af///wEAAAB/\nAAAA/wAAAHn///8B////Af///wH///8BAAAAgQAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8A\nAACZ////Af///wH///8B////Af///wEAAACLAAAA/wAAAP8AAAD/AAAA/wAAAA////8B////AQAA\nABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAh////\nAf///wH///8B////Af///wEAAABDAAAA/wAAAMf///8B////Af///wH///8BAAAAEwAAAO0AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAPcAAABZ////Af///wH///8B////Af///wEAAABvAAAA/wAAAP8A\nAAD/AAAA8////wH///8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAAAh////Af///wH///8B////Af///wEAAAALAAAA/QAAAP0AAAAX////\nAf///wH///8B////AQAAAEUAAAD1AAAA/wAAAP8AAAD/AAAAxQAAACf///8B////Af///wH///8B\n////Af///wEAAAA/AAAAwQAAAMEAAADBAAAAo////wH///8B////AQAAABEAAADBAAAAwQAAAMEA\nAADBAAAAwQAAAMEAAADBAAAAwQAAAMEAAADBAAAAwQAAAMEAAAAZ////Af///wH///8B////Af//\n/wH///8BAAAAzQAAAP8AAABh////Af///wH///8B////Af///wEAAAAhAAAAcQAAAGcAAAAn////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH/\n//8B////Af///wH///8B////Af///wH///8BAAAAkQAAAP8AAACv////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8BAAAAVwAAAP8A\nAAD1AAAACf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wH///8BAAAAGwAAAP8AAAD/AAAAS////wH///8B////Af///wH///8B////Af///wH/\n//8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////AQAAAM0AAADpAAAAh////wEAAAAA\nAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAA\nAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA\n//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD/\n/wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//\nAAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8A\nAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAA\nAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8=\n"""

if __name__ == "__main__":
    try:
        main_gui()
    except Exception as x:
        sys.exit("PASTA GUI is exiting because of an error:\n%s " % str(x))


