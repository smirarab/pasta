"""Main script for SATe GUI on Windows/Mac/Linux
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

import os
import platform
import subprocess
import tempfile
import sys
import time
import wx
from sate import PROGRAM_AUTHOR
from sate import PROGRAM_INSTITUTE
from sate import PROGRAM_DESCRIPTION
from sate import PROGRAM_LICENSE
from sate import PROGRAM_NAME
from sate import PROGRAM_VERSION
from sate import PROGRAM_WEBSITE
from sate import PROGRAM_YEAR
from sate import GLOBAL_DEBUG
from sate import DEFAULT_MAX_MB
from ConfigParser import RawConfigParser
from sate import sate_is_frozen
from sate import sate_home_dir
from sate.configure import get_invoke_run_sate_command
from sate.tools import AlignerClasses
from sate.tools import MergerClasses
from sate.tools import TreeEstimatorClasses
from sate.tools import get_aligner_classes, get_merger_classes, get_tree_estimator_classes
from sate import filemgr

WELCOME_MESSAGE = "%s %s, %s\n\n"% (PROGRAM_NAME, PROGRAM_VERSION, PROGRAM_YEAR)
GRID_VGAP = 8
GRID_HGAP = 8


class SateFrame(wx.Frame):

    def __init__(self, size):
        wx.Frame.__init__(self, None, -1, "SATe - Simultaneous Alignment and Tree Estimation", size=(640,480), style=wx.DEFAULT_FRAME_STYLE)
        self.SetBackgroundColour(wx.LIGHT_GREY)
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetStatusText("SATe Ready!")
        if wx.Platform == '__WXMSW__' or wx.Platform == '__WXMAC__':
            import base64
            import cStringIO
            icon = wx.EmptyIcon()
            icon.CopyFromBitmap(wx.BitmapFromImage(wx.ImageFromStream(cStringIO.StringIO(base64.b64decode(ICO_STR)))))
            self.SetIcon(icon)

        self.ctrls = []
        sizer_all = wx.BoxSizer(wx.VERTICAL)

        #self.fontsize = self.GetFont().GetPointSize()
        #self.title = wx.StaticText(self, -1, "SATe - Simultaneous Alignment and Tree Estimation")
        #self.title.SetFont(wx.Font(self.fontsize+10, wx.SWISS, wx.NORMAL, wx.NORMAL))
        #self.sizer_title = wx.BoxSizer(wx.VERTICAL)
        #self.sizer_title.Add(self.title, 0, wx.ALL, 15)
        #sizer_all.Add(self.sizer_title, 0, wx.ALIGN_CENTER_HORIZONTAL)

        self.sizer_tool_settings = self._create_tools_sizer()
        self.sizer_data = self._create_data_sizer()
        self.sizer_sate_settings = self._create_sate_settings_sizer()
        self.sizer_job_settings = self._create_job_settings_sizer()

        sizer1 = wx.BoxSizer(wx.VERTICAL)

        sizer1.Add(self.sizer_tool_settings, 0, wx.EXPAND|wx.BOTTOM|wx.RIGHT, 5)
        sizer1.Add(self.sizer_data, 0, wx.EXPAND|wx.TOP|wx.RIGHT, 5)

        self.sizer_settings = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_settings.Add(sizer1, 0, wx.EXPAND|wx.ALL, 0)

        sizer2 = wx.BoxSizer(wx.VERTICAL)
        sizer2.Add(self.sizer_job_settings, 0, wx.EXPAND|wx.ALL, 0)
        sizer2.Add(self.sizer_sate_settings, 0, wx.EXPAND|wx.ALL, 0)
        self.sizer_settings.Add(sizer2, 0, wx.EXPAND|wx.ALL, 0)

        sizer_all.Add(self.sizer_settings, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 10)

        self.button = wx.Button(self, label="Start")
        self.log = wx.TextCtrl(self, -1, '', size=(200,120),style=wx.TE_MULTILINE|wx.TE_READONLY|wx.TE_RICH2)
        self.log.AppendText(WELCOME_MESSAGE)
        self.log.AppendText('Running Log (%s %s)\n\n' % (time.strftime("%Y-%m-%d %H:%M:%S"), time.tzname[0]))
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

    def _create_job_settings_sizer(self):
        staticboxsizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "Job Settings"), wx.VERTICAL)
        sizer = wx.GridBagSizer(GRID_VGAP, GRID_HGAP)
        cr = 0
        sizer.Add(wx.StaticText(self, -1, "Job Name"),(cr,0), flag=wx.ALIGN_LEFT )
        self.txt_jobname = wx.TextCtrl(self,-1,'satejob')
        sizer.Add(self.txt_jobname, (cr,1), flag=wx.EXPAND)
        cr += 1
        self.outputdir_btn = wx.Button(self, label="Output Dir." )
        sizer.Add(self.outputdir_btn,(cr,0), flag=wx.ALIGN_LEFT )
        self.txt_outputdir = wx.TextCtrl(self, -1, '', size=(250,9))
        sizer.Add(self.txt_outputdir, (cr,1), flag=wx.EXPAND)
        cr += 1
        sizer.Add(wx.StaticText(self, -1, "CPU(s) Available"), (cr,0), flag=wx.ALIGN_LEFT )
        self.cb_ncpu = wx.ComboBox(self, -1, "1", choices=map(str, range(1,9)), style=wx.CB_READONLY)
        sizer.Add(self.cb_ncpu, (cr,1), flag=wx.EXPAND)
        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Maximum MB"), (cr,0), flag=wx.ALIGN_LEFT )
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

    def _set_custom_sate_settings(self, event):
        self.cb_sate_presets.SetValue("(custom)")

    def OnSatePresets(self, event):
        #wx.MessageBox('Preset Selected (Binding 1)', 'Info')
        preset_selection = self.cb_sate_presets.GetValue()
        if preset_selection != "(custom)":
            self.cb_decomp.SetValue("Longest")
            self.blindmode.SetValue(1)
            self.cb_apply_stop_rule.Enable()
            self.rb_stop2.SetValue(1)
            self.cb_stop1.Disable()
            self.cb_stop2.Enable()
            self.cb_stop2.SetValue("10")
            self.cb_apply_stop_rule.SetValue("After Launch")
            if preset_selection == "SATe-II-simple":
                self.cb_tree_and_alignment.SetValue("Final")
            elif preset_selection == "SATe-II-ML":
                self.cb_tree_and_alignment.SetValue("Best")
            elif preset_selection == "SATe-II-fast":
                self.cb_stop2.SetValue("1")
                self.cb_apply_stop_rule.SetValue("After Blind Mode")
                self.cb_tree_and_alignment.SetValue("Best")
        else:
            pass

    def _create_tools_sizer(self):
        from sate.configure import get_configuration
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
        self.fasttree_prot_models = ["JTT+G20", "JTT+CAT"]

        if GLOBAL_DEBUG:
            defaults = {"Aligner":"PADALIGNER", "Merger":"PADALIGNER", "TreeEstimator":"RANDTREE"}
        else:
            defaults = {"Aligner":"MAFFT", "Merger":"OPAL", "TreeEstimator":"RAXML"}
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
        self.datatype = wx.ComboBox(self, -1, "Nucleotide", (-1, -1), (-1, -1), ["Nucleotide", "Protein"], wx.CB_READONLY)
        self.seq_btn = wx.Button(self, label="Sequence file ..." )
        self.tree_btn = wx.Button(self, label="Tree file (optional) ..." )

        self.txt_seqfn = wx.TextCtrl(self,-1)
        self.txt_treefn = wx.TextCtrl(self,-1)

        self.cb_multilocus = wx.CheckBox(self, -1, "Multi-Locus Data")
        #self.cb_multilocus.Disable()

        sizer.AddMany([ (self.seq_btn, 0, wx.LEFT|wx.EXPAND),
                        (self.txt_seqfn, 0),
                        (wx.StaticText(self, -1, ""), 0, wx.EXPAND),
                        (self.cb_multilocus, 1, wx.EXPAND),
                        (wx.StaticText(self, -1, "Data Type"), 0, wx.ALIGN_CENTER),
                        (self.datatype, 0),
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

    def _create_sate_settings_sizer(self):
        staticboxsizer = wx.StaticBoxSizer(wx.StaticBox(self, -1, "SATe Settings"), wx.VERTICAL)
        sizer = wx.GridBagSizer(GRID_VGAP, GRID_HGAP)

        preset_choices = ["SATe-II-fast", "SATe-II-ML", "SATe-II-simple", "(Custom)",]
        self.cb_sate_presets = wx.ComboBox(self,
                -1,
                "SATe-II-fast",
                choices=preset_choices,
                style=wx.CB_READONLY)

        tree_and_alignment_choices = ["Best", "Final"]
        self.cb_tree_and_alignment = wx.ComboBox(self,
                -1,
                tree_and_alignment_choices[0],
                choices=tree_and_alignment_choices,
                style=wx.CB_READONLY)

        timelimit_list = map(str, [i/100.0 for i in range(1,10)] + [i/10.0 for i in range(1,10)] + range(1,73))
        iterlimit_list = map(str, [1, 5, 10, 20, 50, 100, 200, 500, 1000])
        self.rb_maxsub1 = wx.RadioButton(self, -1, "Fraction", name='frac', style=wx.RB_GROUP)
        self.rb_maxsub2 = wx.RadioButton(self, -1, "Size", name='size')
        self.cb_maxsub1 = wx.ComboBox(self, -1, "20", choices=map(str, range(1,51)), style=wx.CB_READONLY)
        self.cb_maxsub2 = wx.ComboBox(self, -1, "200", choices=map(str, range(1,201)), style=wx.CB_READONLY)

        self.ctrls.extend([self.rb_maxsub1,
                           self.cb_maxsub1,
                           self.rb_maxsub2,
                           self.cb_maxsub2
                           ])

        self.rb_stop1 = wx.RadioButton(self, -1, "Time Limit (h)", name="timelimit", style=wx.RB_GROUP)
        self.rb_stop2 = wx.RadioButton(self, -1, "Iteration Limit", name="iterlimit")
        self.cb_stop1 = wx.ComboBox(self, -1, "24", choices=timelimit_list, style=wx.CB_READONLY)
        self.cb_stop2 = wx.ComboBox(self, -1, "8",choices=iterlimit_list, style=wx.CB_READONLY)
        self.blindmode = wx.CheckBox(self, -1, "Blind Mode Enabled")

        apply_stop_rule_choices = ["After Launch", "After Blind Mode"]
        self.cb_apply_stop_rule = wx.ComboBox(self,
                -1,
                apply_stop_rule_choices[0],
                choices=apply_stop_rule_choices,
                style=wx.CB_READONLY)

        self.blindmode.Value = True
        if not self.blindmode.Value:
            self.cb_apply_stop_rule.SetValue("After Launch")
            self.cb_apply_stop_rule.Disable()

        self.ctrls.extend([self.rb_stop1,
                           self.cb_stop1,
                           self.rb_stop2,
                           self.cb_stop2,
                           self.blindmode,
                           self.cb_apply_stop_rule,
                           ])

        strategy_list = ['Centroid', 'Longest']
        self.cb_decomp = wx.ComboBox(self, -1, "Longest", choices=strategy_list, style=wx.CB_READONLY)

        self.ctrls.append(self.cb_decomp)

        cr = 0

        sizer.Add(wx.StaticText(self, -1, "Quick Set"), (cr, 0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.cb_sate_presets, (cr, 1), flag=wx.EXPAND)

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Max. Subproblem"), (cr,0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.rb_maxsub1, (cr,1), flag=wx.ALIGN_LEFT)
        sizer.Add(self.cb_maxsub1, (cr,2), flag=wx.EXPAND)

        cr += 1
        sizer.Add(self.rb_maxsub2, (cr,1), flag=wx.ALIGN_LEFT)
        sizer.Add(self.cb_maxsub2, (cr,2), flag=wx.EXPAND)

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Decomposition"), (cr,0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.cb_decomp, (cr,1), flag=wx.EXPAND)

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Apply Stop Rule"), (cr,0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.cb_apply_stop_rule, (cr,1), flag=wx.EXPAND)

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Stopping Rule"), (cr,0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.blindmode, (cr,1), flag=wx.EXPAND)

        cr += 1
        sizer.Add(self.rb_stop1, (cr,1), flag=wx.ALIGN_LEFT)
        sizer.Add(self.cb_stop1, (cr,2), flag=wx.EXPAND)

        cr += 1
        sizer.Add(self.rb_stop2, (cr,1), flag=wx.ALIGN_LEFT)
        sizer.Add(self.cb_stop2, (cr,2), flag=wx.EXPAND)

        cr += 1
        sizer.Add(wx.StaticText(self, -1, "Return"), (cr, 0), flag=wx.ALIGN_LEFT )
        sizer.Add(self.cb_tree_and_alignment, (cr, 1), flag=wx.EXPAND)

        self.cb_maxsub1.Disable()
        self.cb_maxsub2.Disable()
        self.rb_maxsub1.Value = True
        self.cb_maxsub1.Enable()

        self.cb_stop1.Disable()
        self.cb_stop2.Enable()
        if self.blindmode.Value:
            self.rb_stop2.Value = True
            self.cb_stop2.Value = "100"

        self.Bind(wx.EVT_COMBOBOX, self.OnSatePresets, self.cb_sate_presets)
        self.OnSatePresets(self.cb_sate_presets)

        self.Bind(wx.EVT_CHECKBOX, self.OnBlindMode, self.blindmode)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnMaxSubproblem, self.rb_maxsub1)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnMaxSubproblem, self.rb_maxsub2)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnStopRule, self.rb_stop1)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnStopRule, self.rb_stop2)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_sate_settings, self.cb_decomp)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_sate_settings, self.cb_apply_stop_rule)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_sate_settings, self.cb_stop1)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_sate_settings, self.cb_stop2)
        self.Bind(wx.EVT_COMBOBOX, self._set_custom_sate_settings, self.cb_tree_and_alignment)

        #cr += 1
        #presets = wx.ComboBox(self, -1, "1", choices=map(str, range(1,9)), style=wx.CB_READONLY)
        #sizer.Add(wx.StaticText(self, -1, "Preset Configuration"), (cr,0), flag=wx.ALIGN_LEFT )
        #sizer.Add(presets, (cr,1), flag=wx.EXPAND)

        staticboxsizer.Add(sizer, 0, wx.ALL, 0)
        return staticboxsizer

    def _create_menu(self):
        self.menuBar = wx.MenuBar()
        self.menuFile = wx.Menu()
        self.menuHelp = wx.Menu()
        self.menuFileSaveLog = self.menuFile.Append(-1, "&Save Log...\tCtrl+S")
        self.menuFileExit = self.menuFile.Append(wx.ID_EXIT, "&Quit SATe\tCtrl+Q")
        self.menuHelpHelp = self.menuHelp.Append( -1, "&Help")
        self.menuHelpAbout = self.menuHelp.Append(wx.ID_ABOUT, "&About SATe")
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
        if self.datatype.Value == "Nucleotide":
            self.cb_tools["model"].Clear()
            if self.cb_tools["treeestimator"].Value.lower() == "raxml":
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
            if not fn[-4:] == '.log':
                fn += ".log"
        else:
            fn += ".log"
        fc = open(fn, 'w')
        fc.write(self.log.GetValue())
        fc.close()

    def OnMaxSubproblem(self, event):
        radio_selected = event.GetEventObject()
        if radio_selected.GetName() == "frac":
            self.cb_maxsub1.Enable()
            self.cb_maxsub2.Disable()
        elif radio_selected.GetName() == "size":
            self.cb_maxsub2.Enable()
            self.cb_maxsub1.Disable()

    def OnBlindMode(self, event):
        self._set_custom_sate_settings(event)
        if self.blindmode.Value:
            self.cb_apply_stop_rule.Enable()
        else:
            self.cb_apply_stop_rule.SetValue("After Launch")
            self.cb_apply_stop_rule.Disable()

    def OnStopRule(self, event):
        self._set_custom_sate_settings(event)
        cb = event.GetEventObject()
        if cb.GetName() == "timelimit":
            self.cb_stop1.Enable()
            self.cb_stop2.Disable()
        elif cb.GetName() == "iterlimit":
            self.cb_stop2.Enable()
            self.cb_stop1.Disable()

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
        info.SetCopyright('Copyright (C) %s' % PROGRAM_YEAR)
        info.SetWebSite((PROGRAM_WEBSITE, '%s Homepage' % PROGRAM_NAME))
        info.SetLicense(PROGRAM_LICENSE)
        info.SetDescription(PROGRAM_DESCRIPTION)
        [info.AddDeveloper(i) for i in PROGRAM_AUTHOR]
        wx.AboutBox(info)

    def OnChooseSeq(self, event):
        if not self.cb_multilocus.Value:
            dialog = wx.FileDialog(None, "Choose sequences...", wildcard = "FASTA files (*.fasta)|*.fasta|FASTA files (*.fas)|*.fas", style=wx.FD_OPEN)
            dialog.ShowModal()
            self.txt_seqfn.SetValue( dialog.GetPath() )
            f = self.txt_seqfn.GetValue()
            if f and not self.txt_outputdir.GetValue():
                self.txt_outputdir.SetValue(os.path.dirname(os.path.abspath(f)))
        else:
            dialog = wx.DirDialog(None, "Choose directory for multiple sequence files", style=wx.FD_OPEN)
            dialog.ShowModal()
            self.txt_seqfn.SetValue( dialog.GetPath() )
            f = self.txt_seqfn.GetValue()
            if f and not self.txt_outputdir.GetValue():
                self.txt_outputdir.SetValue(os.path.abspath(f))

    def OnChooseTree(self, event):
        dialog = wx.FileDialog(None, "Choose tree...", wildcard = "Tree files (*.tree)|*.tree|Tree files (*.tre)|*.tre|Tree files (*.phy)|*.phy", style=wx.FD_OPEN)
        dialog.ShowModal()
        self.txt_treefn.SetValue( dialog.GetPath() )

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
        self.log.AppendText('Job %s is finished.\n' % self.txt_jobname.GetValue())
        self._remove_config_file()
        self._ReactivateOptions()
        self.statusbar.SetStatusText("SATe Ready!")
        self.button.SetLabel('Start')

    def OnButton(self, event):
        if self.button.GetLabel() == "Start":
            self._OnStart()
        elif self.button.GetLabel() == 'Stop':
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
            cfg_success = self._create_config_file()
            if not cfg_success:
                return
            #command = [filemgr.quoted_file_path(x) for x in get_invoke_run_sate_command()]
            command = get_invoke_run_sate_command()
            input_filename = self.txt_seqfn.GetValue()
            if not input_filename or not os.path.isfile(input_filename) and (not self.cb_multilocus.Value):
                wx.MessageBox("Sequence file name is REQUIRED by SATe!", "WARNING", wx.OK|wx.ICON_WARNING)
                self._remove_config_file()
                return
            treefilename = self.txt_treefn.GetValue()
            jobname = self.txt_jobname.GetValue()
            if not jobname:
                wx.MessageBox("Job name cannot be empty, it is REQUIRED by SATe!", "WARNING", wx.OK|wx.ICON_WARNING)
                self._remove_config_file()
                return

            command.extend(['-i', filemgr.quoted_file_path(input_filename)])
            if treefilename and os.path.isfile(treefilename):
                command.extend(['-t', filemgr.quoted_file_path(treefilename)])
            command.extend(['-j', filemgr.quoted_file_path(jobname) ])
            if self.datatype.Value == "Nucleotide":
                dt = "dna"
            else:
                dt = "protein"
            command.extend(['-d', dt])
            command.extend(['%s' % filemgr.quoted_file_path(self.process_cfg_file)])
            self.process = wx.Process(self)
            self.process.Redirect()
            self.pid = wx.Execute( ' '.join(command), wx.EXEC_ASYNC, self.process)
            self.button.SetLabel("Stop")
            self.statusbar.SetStatusText("SATe Running!")
            self._FreezeOptions()

        else:
            self.log.AppendText('Job %s is still running!\n' % self.txt_jobname.GetValue())

    def _OnStop(self):
        if self.process is not None:
            self.log.AppendText('Job %s is terminated early.\n' % self.txt_jobname.GetValue())
            self.process.Kill(self.pid, wx.SIGKILL)
            self._remove_config_file()
            self._ReactivateOptions()
            self.button.SetLabel('Start')
            self.statusbar.SetStatusText("SATe Ready!")
        else:
            self.log.AppendText('No active SATe jobs to terminate!\n')

    def _create_config_file(self):
        from sate.configure import get_configuration
        cfg = get_configuration()

        #if self.txt_resultdir.Value:
        #    basefilename = os.path.basename(self.txt_seqfn.GetValue())
        #    jobname = self.txt_jobname.GetValue()
        #    resultdir = self.txt_resultdir.Value
        #    cfg.commandline.output = os.path.join(resultdir, basefilename+"_%s.aln" % jobname )
        #    cfg.commandline.result = os.path.join(resultdir, basefilename+"_%s.tre" % jobname )

        cfg.sate.aligner = self.cb_tools["aligner"].Value
        cfg.sate.merger = self.cb_tools["merger"].Value
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
            elif model_desc == "JTT":
                cfg.fasttree.model = "-gtr"
            else:
                raise Exception("Unrecognized model: %s" % model_desc)
        cfg.sate.break_strategy = self.cb_decomp.Value
        cfg.sate.start_tree_search_from_current = True
        cfg.commandline.keeptemp = True
        cfg.commandline.keepalignmenttemps = True

        if self.rb_maxsub1.Value:
            cfg.sate.max_subproblem_frac = float(self.cb_maxsub1.Value)/100.0
        elif self.rb_maxsub2.Value:
            cfg.sate.max_subproblem_size = self.cb_maxsub2.Value

        if self.cb_multilocus.Value:
            cfg.commandline.multilocus = True

        if self.blindmode.Value:
            cfg.sate.move_to_blind_on_worse_score = True
            if self.cb_apply_stop_rule.GetValue() == "After Blind Mode":
                if self.rb_stop1.Value:
                    cfg.sate.after_blind_time_without_imp_limit = float(self.cb_stop1.Value)*3600
                elif self.rb_stop2.Value:
                    cfg.sate.after_blind_iter_without_imp_limit = float(self.cb_stop2.Value)
            else:
                if self.rb_stop1.Value:
                    cfg.sate.time_limit = float(self.cb_stop1.Value)*3600
                elif self.rb_stop2.Value:
                    cfg.sate.iter_limit = self.cb_stop2.Value
        else:
            if self.rb_stop1.Value:
                cfg.sate.time_limit = float(self.cb_stop1.Value)*3600
            elif self.rb_stop2.Value:
                cfg.sate.iter_limit = self.cb_stop2.Value

        cfg.sate.return_final_tree_and_alignment = self.cb_tree_and_alignment.GetValue() == "Final"
        cfg.sate.output_directory = self.txt_outputdir.GetValue()
        cfg.sate.num_cpus = self.cb_ncpu.Value
        max_mb = self.txt_maxmb.GetValue()
        if not self.validate_max_mb(max_mb):
            return False
        cfg.sate.max_mem_mb = max_mb

        # this creates a file that cannot be deleted while the Python
        # process is running (under the mess that is called 'Windows')
        #tf, self.process_cfg_file = tempfile.mkstemp(dir=sate_home_dir(),
        #        suffix='_internal.cfg')

        tf = tempfile.NamedTemporaryFile(suffix='_internal.cfg', dir=sate_home_dir())
        self.process_cfg_file = tf.name
        tf.close()
        cfg.save_to_filepath(self.process_cfg_file)
        return True

    def _remove_config_file(self):
        if "SATE_GUIDEVMODE" in os.environ:
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

class SateApp(wx.PySimpleApp):
    def OnInit(self):
        self.frame = SateFrame(size=wx.Display().GetClientArea())
        self.frame.Show(True)
        self.SetTopWindow(self.frame)
        return True

def main_gui():
    app = SateApp()
    app.MainLoop()

ICO_STR = """AAABAAMAEBAAAAAAIABoBAAANgAAACAgAAAAACAAqBAAAJ6EAAAwMAAAAAAgAKglAABGFQAAKAAA\nABAAAAAgAAAAAQAgAAAAAABABAAAAAAAAAAAAAAAAAAAAAAAAAAAAGsAAADvAAAAqQAAAEUAAACX\n////Af///wEAAAC3AAAAJQAAAGsAAABv////Af///wEAAACHAAAA8QAAAGkAAAD1AAAA/wAAAP8A\nAABpAAAA4f///wEAAAALAAAA/wAAABkAAACPAAAAk////wEAAAAVAAAA+wAAAP8AAADXAAAA8wAA\nALsAAAD7AAAAgQAAAPsAAAAlAAAAQwAAAPkAAAADAAAAjwAAAJP///8BAAAAUQAAAO0AAABfAAAA\ntwAAAGv///8BAAAAsQAAAHsAAAD/AAAA/wAAAP8AAADd////AQAAAI8AAACT////AQAAAHUAAACl\n////AQAAAB3///8B////AQAAAKcAAAB7AAAA6QAAAP8AAAD/AAAAv////wEAAACPAAAAk////wEA\nAACHAAAAsQAAAFMAAABT////AQAAABcAAADlAAAAcwAAAM0AAACvAAAAwQAAAKP///8BAAAAjwAA\nAJP///8BAAAAjwAAAP8AAAD/AAAA/wAAAB0AAADfAAAA/wAAAFsAAACvAAAAawAAAJMAAACH////\nAQAAAI8AAACT////AQAAAIsAAADnAAAAzwAAAPkAAACZAAAA/wAAAP0AAAAhAAAAkwAAAIUAAACv\nAAAAaf///wEAAACPAAAAk////wEAAAB7AAAAmQAAACMAAADrAAAA2QAAAP8AAACF////AQAAAHUA\nAAChAAAAyQAAAE3///8BAAAAjwAAAJP///8BAAAAWwAAANUAAABjAAAAzQAAAPMAAABv////Af//\n/wEAAABZAAAAvQAAAOUAAAAv////AQAAAI8AAACT////AQAAACMAAAD/AAAA/wAAAJUAAAD9AAAA\nE////wH///8BAAAAOwAAANsAAAD7AAAAEf///wEAAACPAAAAk////wH///8BAAAAtwAAAP0AAAAz\nAAAA9QAAACsAAAA3AAAAKQAAAB8AAAD9AAAA8wAAAAMAAAA3AAAApwAAAKkAAAA3AAAAAwAAAAsA\nAAAp////AQAAANkAAADvAAAA/QAAADMAAAAFAAAA+wAAANcAAAAJAAAA/wAAAP8AAAD/AAAA/wAA\nAAsAAAAFAAAAXf///wEAAACdAAAA/wAAAP8AAAAz////AQAAAOMAAAC5AAAACQAAAP8AAAD/AAAA\n/wAAAP8AAAAL////AQAAAKf///8BAAAAJQAAAMUAAACPAAAAC////wEAAAB3AAAAXwAAAAUAAACV\nAAAAlQAAAJUAAACVAAAAB////wEAAACZAAAAIf///wH///8B////Af///wH///8B////Af///wH/\n//8B////Af///wH///8B////Af///wH///8BAAAAXQAAAGsAAP//AAD//wAA//8AAP//AAD//wAA\n//8AAP//AAD//wAA//8AAP//AAD//wAA//8AAP//AAD//wAA//8AAP//KAAAACAAAABAAAAAAQAg\nAAAAAACAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAkAAABTAAAAyQAAAPkAAAC7AAAAMwAAAAcAAABb\nAAAAgQAAAEX///8B////Af///wH///8BAAAAbQAAAIEAAAA3////AQAAAB0AAABzAAAAdQAAAB//\n//8B////Af///wH///8BAAAAFwAAAJUAAAD3AAAA0QAAAFUAAAAJAAAAZwAAAOcAAAD/AAAA/wAA\nAP0AAAC1AAAABQAAAK0AAAD/AAAAmf///wH///8B////Af///wEAAADrAAAA/wAAAF3///8BAAAA\nOwAAAOUAAADnAAAAP////wH///8B////Af///wEAAAB5AAAA9wAAAP8AAAD/AAAA5wAAAF8AAADp\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAxAAAAkQAAAP8AAAC3////Af///wH///8BAAAACwAAAP8A\nAAD/AAAAPf///wEAAAA7AAAA5QAAAOcAAAA/////Af///wH///8BAAAADQAAAO0AAAD/AAAA/wAA\nAP8AAAD/AAAArwAAAOkAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAHEAAAB3AAAA/wAAAM////8B////\nAf///wEAAAAjAAAA/wAAAP8AAAAj////AQAAADsAAADlAAAA5wAAAD////8B////Af///wEAAABF\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAACvAAAA6QAAAP8AAAD/AAAA/wAAAP8AAAD/AAAApwAAAFkA\nAAD/AAAA7////wH///8B////AQAAAEEAAAD/AAAA+wAAAAf///8BAAAAOwAAAOUAAADnAAAAP///\n/wH///8B////AQAAAI0AAAD/AAAA+wAAAKsAAAC3AAAA9QAAAK8AAADpAAAA+QAAAIUAAABpAAAA\n8QAAAP8AAAC5AAAASwAAAP8AAAD9AAAASQAAAEcAAABHAAAAgwAAAP8AAADp////Af///wEAAAA7\nAAAA5QAAAOcAAAA/////Af///wEAAAAHAAAArQAAAP8AAAC/AAAACQAAABMAAACPAAAAqwAAANUA\nAABX////Af///wEAAAB7AAAA/wAAAMUAAAA3AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAMn///8B////AQAAADsAAADlAAAA5wAAAD////8B////AQAAABsAAADFAAAA/QAAAGP///8B////\nAQAAABMAAABXAAAAeQAAAAv///8B////AQAAAFMAAAD5AAAAywAAACcAAAD7AAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAArf///wH///8BAAAAOwAAAOUAAADnAAAAP////wH///8BAAAAJQAAAM8A\nAADvAAAARf///wH///8B////AQAAAAkAAAAF////Af///wH///8BAAAATQAAAPcAAADPAAAAJwAA\nAOEAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACN////Af///wEAAAA7AAAA5QAAAOcAAAA/////\nAf///wEAAAAxAAAA2wAAAOMAAAA5////Af///wH///8B////Af///wH///8B////Af///wEAAABd\nAAAA/QAAANEAAAAnAAAAxwAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAHP///8B////AQAAADsA\nAADlAAAA5wAAAD////8B////AQAAADUAAADfAAAA8wAAALkAAACnAAAApwAAAKcAAACl////Af//\n/wH///8BAAAAAwAAAJ8AAAD/AAAAywAAACEAAACnAAAA/wAAAPEAAADbAAAA3wAAAPkAAAD9AAAA\nVf///wH///8BAAAAOwAAAOUAAADnAAAAP////wH///8BAAAAOQAAAOMAAAD/AAAA/wAAAP8AAAD/\nAAAA/wAAAP////8B////AQAAAAMAAABTAAAA9QAAAP8AAADFAAAAGwAAAI0AAAD/AAAAuwAAADMA\nAABTAAAA3QAAAPEAAABJ////Af///wEAAAA7AAAA5QAAAOcAAAA/////Af///wEAAAA5AAAA4wAA\nAP8AAAD/AAAA/wAAAP8AAAD/AAAA/////wEAAAAHAAAAjwAAAPUAAAD/AAAA/wAAALkAAAARAAAA\nawAAAP8AAAC5AAAADwAAADkAAADhAAAA4QAAADf///8B////AQAAADsAAADlAAAA5wAAAD////8B\n////AQAAADcAAADhAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD3AAAAAwAAAGkAAAD5AAAA/wAAAP8A\nAAD/AAAAoQAAAAMAAABVAAAA+wAAAMUAAAAbAAAARQAAAO0AAADVAAAAK////wH///8BAAAAOwAA\nAOUAAADnAAAAP////wH///8BAAAAMwAAAN0AAADxAAAAqwAAAJUAAACrAAAA/wAAAPEAAAAvAAAA\n1wAAAP8AAAD/AAAA/wAAAP8AAABn////AQAAAEUAAADtAAAA1QAAACsAAABVAAAA+wAAAMMAAAAb\n////Af///wEAAAA7AAAA5QAAAOcAAAA/////Af///wEAAAArAAAA0wAAAOcAAAA/////AQAAADsA\nAAD/AAAA3wAAAGUAAAD7AAAA/wAAAP8AAAD/AAAA9wAAAB3///8BAAAANwAAAOEAAADhAAAANwAA\nAGkAAAD/AAAAtwAAAA3///8B////AQAAADsAAADlAAAA5wAAAD////8B////AQAAACEAAADLAAAA\n9QAAAEv///8BAAAATwAAAP8AAADPAAAAnwAAAP8AAAD/AAAA/wAAAPUAAAB/////Af///wEAAAAn\nAAAA0QAAAPEAAABHAAAAhwAAAP8AAACj////Af///wH///8BAAAAOwAAAOUAAADnAAAAP////wH/\n//8BAAAADwAAALkAAAD/AAAAfwAAAAMAAACFAAAA/wAAAKsAAADLAAAA/wAAAP8AAAD/AAAAiQAA\nABP///8B////AQAAABsAAADDAAAA+wAAAFMAAACfAAAA/wAAAIv///8B////Af///wEAAAA7AAAA\n5QAAAOcAAAA/////Af///wEAAAADAAAAowAAAP8AAADXAAAAPQAAAMcAAAD/AAAAhwAAAOMAAAD/\nAAAA4wAAAFv///8B////Af///wH///8BAAAACQAAALMAAAD/AAAAbQAAAL0AAAD/AAAAa////wH/\n//8B////AQAAADsAAADlAAAA5wAAAD////8B////Af///wEAAABlAAAA/wAAAP8AAAD/AAAA/wAA\nAPcAAABTAAAA7QAAAP8AAAB3AAAAB////wH///8B////Af///wH///8BAAAAowAAAP8AAACHAAAA\n1QAAAP8AAABR////Af///wH///8BAAAAOwAAAOUAAADnAAAAP////wH///8B////AQAAACUAAAD/\nAAAA/wAAAP8AAAD/AAAA2wAAADEAAAD7AAAA/wAAAC////8B////Af///wH///8B////Af///wEA\nAACDAAAA/wAAAKUAAADzAAAA/wAAAC////8B////Af///wEAAAA7AAAA5QAAAOcAAAA/////Af//\n/wH///8B////AQAAALEAAAD/AAAA/wAAAP8AAACZAAAABwAAAPsAAAD/AAAAGf///wH///8B////\nAf///wH///8B////AQAAAGkAAAD/AAAAywAAAP8AAAD/AAAAFf///wH///8B////AQAAADsAAADl\nAAAA5wAAAD////8B////Af///wH///8BAAAARwAAAOUAAAD/AAAA8wAAAC3///8BAAAA9QAAAP8A\nAAAx////Af///wEAAAAdAAAAPf///wH///8BAAAASQAAAP8AAAD7AAAA/wAAAPX///8B////Af//\n/wH///8BAAAAOwAAAOUAAADnAAAAP////wH///8B////Af///wH///8BAAAAKQAAAHMAAAAx////\nAf///wEAAADlAAAA/wAAAHUAAAADAAAAGwAAAKUAAABn////Af///wEAAAAvAAAA/wAAAP8AAAD/\nAAAA2////wEAAAAHAAAAawAAAGsAAACNAAAA7wAAAPEAAACPAAAAawAAAGsAAAAJ////Af///wH/\n//8B////Af///wH///8B////AQAAAMcAAAD/AAAA9QAAAMkAAAD1AAAA/wAAAGf///8B////AQAA\nAA8AAAD/AAAA/wAAAP8AAAC7////AQAAAA8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAABf///8B////AQAAAAkAAABDAAAAD////wH///8BAAAAowAAAP8AAAD/AAAA/wAAAP8AAAD/\nAAAAZ////wH///8B////AQAAAPMAAAD/AAAA/wAAAKH///8BAAAADwAAAP8AAAD/AAAA/wAAAP8A\nAAD/AAAA/wAAAP8AAAD/AAAAF////wH///8BAAAACwAAAN8AAABB////Af///wEAAABlAAAA+wAA\nAP8AAAD/AAAA/wAAAP8AAABn////Af///wH///8BAAAA0wAAAP8AAAD/AAAAf////wEAAAAPAAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAX////Af///wH///8BAAAAwwAAAIP///8B\n////AQAAADcAAADfAAAA/wAAAP8AAAD/AAAA/wAAAGf///8B////Af///wEAAAC5AAAA/wAAAP8A\nAABl////AQAAAA8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAABf///8B////Af//\n/wEAAACPAAAAyf///wH///8BAAAABQAAAH0AAAD9AAAA/wAAAPkAAADFAAAAKf///wH///8B////\nAQAAAI0AAADrAAAA5QAAAEn///8BAAAADwAAAOsAAADrAAAA6wAAAOsAAADrAAAA6wAAAOsAAADr\nAAAAFf///wH///8B////AQAAAFUAAAD5AAAAIf///wH///8BAAAAEQAAAHsAAACbAAAAYQAAAB//\n//8B////Af///wH///8BAAAAJQAAAEEAAAA9AAAAE////wEAAAAFAAAAQQAAAEEAAABBAAAAQQAA\nAEEAAABBAAAAQQAAAEEAAAAH////Af///wH///8BAAAANwAAAOEAAABj////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wEAAAAXAAAAwQAAAK8A\nAAAL////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAQAAAAMAAACXAAAAxwAAACcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACgAAAAwAAAAYAAAAAEAIAAAAAAAgCUAAAAA\nAAAAAAAAAAAAAAAAAAD///8B////AQAAAC0AAACTAAAA4QAAAP8AAADNAAAAOf///wH///8BAAAA\nCwAAAEEAAABBAAAAQQAAABH///8B////Af///wH///8B////Af///wEAAAAvAAAAQQAAAEEAAAAt\n////Af///wH///8BAAAAKwAAAEEAAABBAAAAL////wH///8B////Af///wH///8B////Af///wH/\n//8BAAAACwAAAJsAAAD3AAAA7QAAAJ8AAAAx////Af///wH///8BAAAARwAAAPUAAAD/AAAA/wAA\nAP8AAAD/AAAA9QAAADH///8BAAAAHQAAAP8AAAD/AAAA/wAAAFP///8B////Af///wH///8B////\nAf///wEAAADPAAAA/wAAAP8AAACf////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B\n////Af///wH///8B////Af///wH///8BAAAAtwAAAP8AAAD/AAAA/wAAAP8AAAD1AAAARf///wEA\nAABZAAAA+wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAMH///8BAAAAAwAAAPsAAAD/AAAA/wAA\nAG////8B////Af///wH///8B////Af///wEAAADrAAAA/wAAAP8AAACD////Af///wH///8BAAAA\nrwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////Af///wEAAABjAAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAA+QAAAEcAAADfAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8A\nAAAv////AQAAAOEAAAD/AAAA/wAAAIn///8B////Af///wH///8B////AQAAAAcAAAD/AAAA/wAA\nAP8AAABl////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////\nAf///wEAAADVAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAIcAAADfAAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAB7////AQAAAMUAAAD/AAAA/wAAAKX///8B////Af///wH/\n//8B////AQAAACEAAAD/AAAA/wAAAP8AAABJ////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf//\n/wH///8B////Af///wH///8B////AQAAADkAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAAIcAAADfAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAC/////AQAAAKcAAAD/\nAAAA/wAAAMH///8B////Af///wH///8B////AQAAAD0AAAD/AAAA/wAAAP8AAAAr////Af///wH/\n//8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////AQAAAH8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAIcAAADfAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAAP8AAADv////AQAAAIsAAAD/AAAA/wAAAN3///8B////Af///wH///8B////AQAAAFkAAAD/\nAAAA/wAAAP8AAAAN////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH/\n//8B////AQAAAMMAAAD/AAAA/wAAAP8AAAD9AAAA5QAAAP8AAAD/AAAA/wAAAIcAAADfAAAA/wAA\nAP8AAAD/AAAA/QAAAP8AAAD/AAAA/wAAAP8AAAD/AAAAHQAAAG0AAAD/AAAA/wAAAPf///8B////\nAf///wH///8B////AQAAAHUAAAD/AAAA/wAAAPH///8B////Af///wH///8BAAAArwAAAP8AAAD/\nAAAAuf///wH///8B////Af///wH///8B////AQAAAPMAAAD/AAAA/wAAAO8AAAAhAAAAAwAAAFMA\nAADVAAAA/wAAAIcAAADfAAAA/wAAAPMAAABlAAAACwAAACcAAADhAAAA/wAAAP8AAAD/AAAAOQAA\nAE8AAAD/AAAA/wAAAP8AAABxAAAAaQAAAGkAAABpAAAAaQAAALkAAAD/AAAA/wAAANP///8B////\nAf///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAHQAAAP8AAAD/\nAAAA/wAAAHX///8B////Af///wEAAAAfAAAA7wAAAIcAAADfAAAA9QAAADX///8B////Af///wEA\nAABTAAAA/wAAAP8AAAD/AAAASwAAADMAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAALf///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////\nAf///wH///8BAAAARwAAAP8AAAD/AAAA/wAAACX///8B////Af///wH///8BAAAAUQAAAIcAAADf\nAAAAVf///wH///8B////Af///wEAAAALAAAA+wAAAP8AAAD/AAAAXwAAABUAAAD/AAAA/wAAAP8A\nAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAJn///8B////Af///wH///8BAAAArwAA\nAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAYQAAAP8AAAD/AAAA5////wH///8B////\nAf///wH///8B////AQAAACkAAACL////Af///wH///8B////Af///wH///8BAAAA6QAAAP8AAAD/\nAAAAZ////wEAAAD3AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAHv/\n//8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAdwAA\nAP8AAAD/AAAAw////wH///8B////Af///wH///8B////Af///wEAAAAJ////Af///wH///8B////\nAf///wH///8BAAAA3QAAAP8AAAD/AAAAb////wEAAADbAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAF////8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH/\n//8B////Af///wH///8BAAAAiwAAAP8AAAD/AAAAr////wH///8B////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8BAAAA8wAAAP8AAAD/AAAAd////wEAAAC9AAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAEH///8B////Af///wH///8B\nAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAmwAAAP8AAAD/AAAAm////wH/\n//8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wEAAAAVAAAA/wAA\nAP8AAAD/AAAAc////wEAAAChAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAACX///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B\nAAAAoQAAAP8AAAD/AAAA/QAAAPkAAAD5AAAA+QAAAPkAAAD5AAAA+QAAAPf///8B////Af///wH/\n//8B////Af///wEAAAA/AAAA/wAAAP8AAAD/AAAAZf///wEAAACDAAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAAf///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAA\nuf///wH///8B////Af///wH///8BAAAApwAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/\nAAAA/wAAAP////8B////Af///wH///8B////AQAAAA0AAADTAAAA/wAAAP8AAAD/AAAAV////wEA\nAABnAAAA/wAAAP8AAAD7AAAAkQAAAJEAAACRAAAAwQAAAP8AAAD/AAAA6f///wH///8B////Af//\n/wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAArQAAAP8AAAD/AAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP////8B////Af///wH///8BAAAABwAAALEAAAD/\nAAAA/wAAAP8AAAD/AAAAS////wEAAABJAAAA/wAAAP8AAAD/AAAACf///wH///8BAAAAgwAAAP8A\nAAD/AAAAzf///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af//\n/wH///8BAAAAqwAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP3///8B////\nAf///wEAAAAlAAAA0wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAAO////wEAAAArAAAA/wAAAP8AAAD/\nAAAAI////wH///8BAAAAnwAAAP8AAAD/AAAAr////wH///8B////Af///wH///8BAAAArwAAAP8A\nAAD/AAAAuf///wH///8B////Af///wH///8BAAAApQAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAPX///8B////AQAAAB8AAADjAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\nFf///wEAAAAPAAAA/wAAAP8AAAD/AAAAP////wH///8BAAAAuQAAAP8AAAD/AAAAk////wH///8B\n////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAnwAAAP8A\nAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAO////8BAAAABwAAANsAAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAADj////Af///wH///8BAAAA8QAAAP8AAAD/AAAAW////wH///8BAAAA\n1QAAAP8AAAD/AAAAdf///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B\n////Af///wH///8BAAAAlwAAAP8AAAD/AAAAwwAAAGEAAABhAAAAYQAAAI8AAAD/AAAA/wAAAOf/\n//8BAAAAYwAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACx////Af///wH///8BAAAA1QAA\nAP8AAAD/AAAAd////wH///8BAAAA8QAAAP8AAAD/AAAAV////wH///8B////Af///wH///8BAAAA\nrwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAhQAAAP8AAAD/AAAAsf///wH///8B\n////AQAAAFUAAAD/AAAA/wAAANf///8BAAAA1wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8A\nAABv////Af///wH///8BAAAAtwAAAP8AAAD/AAAAkf///wEAAAALAAAA/wAAAP8AAAD/AAAAO///\n/wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAA\nbwAAAP8AAAD/AAAAyf///wH///8B////AQAAAGEAAAD/AAAA/wAAAMMAAAAtAAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAOsAAAAL////Af///wH///8BAAAAmQAAAP8AAAD/AAAArf///wEA\nAAAnAAAA/wAAAP8AAAD/AAAAHf///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf//\n/wH///8B////Af///wH///8BAAAAWQAAAP8AAAD/AAAA6////wH///8B////AQAAAH8AAAD/AAAA\n/wAAAK8AAABfAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAHf///8B////Af///wH///8B\nAAAAfQAAAP8AAAD/AAAAyf///wEAAABBAAAA/wAAAP8AAAD9AAAAA////wH///8B////Af///wH/\n//8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8BAAAAOQAAAP8AAAD/AAAA/wAA\nACX///8B////AQAAALEAAAD/AAAA/wAAAI8AAACRAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA\nnwAAAAP///8B////Af///wH///8BAAAAXwAAAP8AAAD/AAAA5f///wEAAABdAAAA/wAAAP8AAADj\n////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH/\n//8BAAAAEQAAAP8AAAD/AAAA/wAAAHf///8BAAAACwAAAO8AAAD/AAAA/wAAAGcAAAC/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAACN////Af///wH///8B////Af///wH///8BAAAAQwAAAP8AAAD/AAAA\n/QAAAAMAAAB3AAAA/wAAAP8AAADF////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/\nAAAAuf///wH///8B////Af///wH///8B////AQAAAOsAAAD/AAAA/wAAAOsAAAA/AAAAkQAAAP8A\nAAD/AAAA/wAAAD8AAADRAAAA/wAAAP8AAAD/AAAA7QAAAEn///8B////Af///wH///8B////Af//\n/wH///8BAAAAJQAAAP8AAAD/AAAA/wAAABsAAACTAAAA/wAAAP8AAACp////Af///wH///8B////\nAf///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////AQAAAK8AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA+wAAAAsAAADdAAAA/wAAAP8AAAD3AAAAO////wH/\n//8B////Af///wH///8B////Af///wH///8BAAAACQAAAP8AAAD/AAAA/wAAADcAAACvAAAA/wAA\nAP8AAACL////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////\nAf///wH///8B////AQAAAGsAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAAxf///wEAAADp\nAAAA/wAAAP8AAACB////Af///wH///8B////Af///wH///8B////Af///wH///8B////AQAAAOsA\nAAD/AAAA/wAAAFMAAADJAAAA/wAAAP8AAABv////Af///wH///8B////Af///wH///8BAAAArwAA\nAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////AQAAAB8AAAD9AAAA/wAAAP8AAAD/AAAA\n/wAAAP8AAAD/AAAAef///wEAAAD3AAAA/wAAAP8AAABT////Af///wH///8B////Af///wH///8B\n////Af///wH///8B////AQAAAM0AAAD/AAAA/wAAAG8AAADlAAAA/wAAAP8AAABR////Af///wH/\n//8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////Af//\n/wEAAACzAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD7AAAAGf///wEAAAD9AAAA/wAAAP8AAAAt////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////AQAAAK8AAAD/AAAA/wAAAI0AAAD9\nAAAA/wAAAP8AAAAz////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH/\n//8B////Af///wH///8B////Af///wEAAAA5AAAA/QAAAP8AAAD/AAAA/wAAAP8AAACf////Af//\n/wEAAAD3AAAA/wAAAP8AAAAj////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAQAAAJMAAAD/AAAA/wAAAMEAAAD/AAAA/wAAAP8AAAAX////Af///wH///8B////Af///wH///8B\nAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B////Af///wH///8BAAAAiQAAAP8A\nAAD/AAAA/wAAAOUAAAAV////Af///wEAAADzAAAA/wAAAP8AAAA3////Af///wH///8B////Af//\n/wEAAAA7////Af///wH///8B////AQAAAHUAAAD/AAAA/wAAAPUAAAD/AAAA/wAAAPn///8B////\nAf///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAAuf///wH///8B////Af///wH///8B\n////Af///wH///8B////AQAAAF0AAACzAAAAnQAAAB3///8B////Af///wEAAADnAAAA/wAAAP8A\nAABr////Af///wH///8B////AQAAAHsAAACZ////Af///wH///8B////AQAAAFkAAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAN3///8B////Af///wH///8B////Af///wH///8BAAAArwAAAP8AAAD/AAAA\nuf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wEAAADPAAAA/wAAAP8AAADPAAAACf///wEAAAAFAAAAbQAAAP8AAACZ////Af///wH/\n//8B////AQAAADsAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAL////8B////AQAAAA0AAAChAAAAoQAA\nAKEAAAChAAAA4QAAAP8AAAD/AAAA5QAAAKEAAAChAAAAoQAAAKEAAAAV////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wEAAAC3AAAA/wAAAP8AAAD/AAAA0wAAAJ0AAADp\nAAAA/wAAAP8AAACZ////Af///wH///8B////AQAAAB0AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAKH/\n//8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAAh////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wEAAACTAAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACZ////Af///wH///8B////AQAAAAUAAAD9\nAAAA/wAAAP8AAAD/AAAA/wAAAIX///8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8A\nAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAh////Af///wH///8B////AQAAACkAAADpAAAAhf//\n/wH///8B////Af///wEAAABlAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACZ////\nAf///wH///8B////Af///wEAAADjAAAA/wAAAP8AAAD/AAAA/wAAAGf///8B////AQAAABcAAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAh////Af///wH/\n//8B////AQAAAAMAAADxAAAA3f///wH///8B////Af///wEAAAAtAAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAACZ////Af///wH///8B////Af///wEAAADHAAAA/wAAAP8AAAD/AAAA\n/wAAAEn///8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/\nAAAA/wAAAP8AAAAh////Af///wH///8B////Af///wEAAAC5AAAA/wAAACv///8B////Af///wH/\n//8BAAAA4QAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAACZ////Af///wH///8B////Af//\n/wEAAACpAAAA/wAAAP8AAAD/AAAA/wAAAC3///8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA\n/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAh////Af///wH///8B////Af///wEAAAB/\nAAAA/wAAAHn///8B////Af///wH///8BAAAAgQAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8A\nAACZ////Af///wH///8B////Af///wEAAACLAAAA/wAAAP8AAAD/AAAA/wAAAA////8B////AQAA\nABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAAh////\nAf///wH///8B////Af///wEAAABDAAAA/wAAAMf///8B////Af///wH///8BAAAAEwAAAO0AAAD/\nAAAA/wAAAP8AAAD/AAAA/wAAAPcAAABZ////Af///wH///8B////Af///wEAAABvAAAA/wAAAP8A\nAAD/AAAA8////wH///8B////AQAAABcAAAD/AAAA/wAAAP8AAAD/AAAA/wAAAP8AAAD/AAAA/wAA\nAP8AAAD/AAAA/wAAAP8AAAAh////Af///wH///8B////Af///wEAAAALAAAA/QAAAP0AAAAX////\nAf///wH///8B////AQAAAEUAAAD1AAAA/wAAAP8AAAD/AAAAxQAAACf///8B////Af///wH///8B\n////Af///wEAAAA/AAAAwQAAAMEAAADBAAAAo////wH///8B////AQAAABEAAADBAAAAwQAAAMEA\nAADBAAAAwQAAAMEAAADBAAAAwQAAAMEAAADBAAAAwQAAAMEAAAAZ////Af///wH///8B////Af//\n/wH///8BAAAAzQAAAP8AAABh////Af///wH///8B////Af///wEAAAAhAAAAcQAAAGcAAAAn////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH/\n//8B////Af///wH///8B////Af///wH///8BAAAAkQAAAP8AAACv////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8BAAAAVwAAAP8A\nAAD1AAAACf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B\n////Af///wH///8BAAAAGwAAAP8AAAD/AAAAS////wH///8B////Af///wH///8B////Af///wH/\n//8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af//\n/wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////Af///wH///8B////\nAf///wH///8B////Af///wH///8B////Af///wH///8B////AQAAAM0AAADpAAAAh////wEAAAAA\nAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAA\nAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA\n//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD/\n/wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//\nAAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8A\nAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8AAAAAAAD//wAA\nAAAAAP//AAAAAAAA//8AAAAAAAD//wAAAAAAAP//AAAAAAAA//8=\n"""

if __name__ == "__main__":
    try:
        main_gui()
    except Exception, x:
        sys.exit("SATe GUI is exiting because of an error:\n%s " % str(x))


