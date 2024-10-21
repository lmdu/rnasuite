import os
import json
import traceback
import multiprocessing

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *
from PySide6.QtSvgWidgets import *

from plot import *
from renv import *
from config import *
from errors import *
from config import *
from tables import *
from workers import *
from dialogs import *
from widgets import *

__all__ = ['RNASuiteMainWindow']

class RNASuiteMainWindow(QMainWindow):
	def __init__(self):
		super().__init__()
		self.setWindowIcon(QIcon('icons/logo.svg'))
		self.setWindowTitle("RNASuite v{}".format(RNASUITE_VERSION))
		self.pyconn, self.rconn = multiprocessing.Pipe()

		self.global_params = {}

		self.create_main_widgets()
		self.create_input_dock()
		#self.create_plot_dock()
		self.create_output_dock()

		self.create_actions()
		self.create_menus()
		self.create_toolbar()
		self.create_statusbar()
		
		self.create_table_manager()
		self.create_r_enviroment()
		self.create_error_prompter()
		self.read_settings()
		self.show()

	def read_settings(self):
		settings = QSettings()
		settings.beginGroup('Window')
		self.resize(settings.value('size', QSize(1000, 700)))
		self.move(settings.value('pos', QPoint(200, 200)))
		settings.endGroup()

	def write_settings(self):
		settings = QSettings()

		if not self.isMaximized():
			settings.beginGroup('Window')
			settings.setValue('size', self.size())
			settings.setValue('pos', self.pos())
			settings.endGroup()

	def create_r_enviroment(self):
		self.r_thread = RNASuiteRMessageProcessor(self)
		self.r_thread.error.connect(self.show_error_message)
		self.r_thread.warning.connect(self.show_warn_message)
		self.r_thread.message.connect(self.show_status_message)
		self.r_thread.socket.connect(self.plot_viewer.connect_to_socket)
		self.r_thread.plot.connect(self.plot_viewer.update_image)
		self.r_thread.results.connect(self.output_list.receive)
		self.r_thread.running.connect(self.wait_spinner.toggle)
		self.r_thread.running.connect(self.run_progress.setVisible)
		self.r_thread.start()

		self.r_process = RNASuiteREnvironment(self.rconn, RNASUITE_PORT, RNASUITE_TOKEN)
		self.r_process.start()

	def create_error_prompter(self):
		self.main_error = RNASuiteErrorPrompter(self)
		self.main_error.error_occurred.connect(self.show_error_message)

	def create_table_manager(self):
		self.table_widgets = RNASuiteTableWidgets(self)
		#self.table_widgets.intab.connect(self.input_tabs.addTab)
		#self.table_widgets.outab.connect(self.output_tabs.addTab)
		#self.table_widgets.title.connect(self.on_change_tab_title)

	@Slot()
	def on_change_tab_title(self, widget, title):
		index = self.output_tabs.indexOf(widget)
		self.output_tabs.setTabText(index, title)

	@Slot()
	def process_r_results(self, rtype, res):
		if rtype == 'degs':
			v = res['degs_versus']
			contrast = "{} vs {}".format(v[-2], v[-1])

			if 'normal_count' in res:
				self.output_list.add_row(
					type = 'Counts',
					name = 'Noramlized Counts',
					plot = 0,
					update = 1,
					data = res['normal_count']
				)

			elif 'degs_list' in res:
				self.output_list.add_row(
					type = 'DEGs',
					name = contrast,
					plot = 0,
					update = 1,
					data = res['degs_list']
				)

	@Slot()
	def on_open_project(self):
		code = """
		x<-c(1,2,3);
		y<-c(4,5,6);
		plot(x, y, type='l')
		"""
		self.pyconn.send({
			'action': 'eval',
			'code': code
		})

	@Slot()
	def on_save_project(self):
		data = self.exp_table.get_data()
		#R.rcall(('convert_pandas_to_dataframe', data), _asis=True)

	@Slot()
	def on_save_project_as(self):
		#R.rcall('test_data', ['groups', 'sex', 'age'])
		pass

	@Slot()
	def on_close_project(self):
		pass

	@Slot()
	def on_export_file(self):
		pass

	@Slot()
	def on_clicked(self):
		text = self.input.text()

		self.pyconn.send({
			'action': 'eval',
			'code': text,
		})

		self.input.setText('')
		#print(R.rcall('ls'))

	def closeEvent(self, event):
		self.write_settings()
		self.plot_viewer.disconnect_to_socket()
		self.pyconn.close()
		event.accept()

	def resizeEvent(self, event):
		size = event.size()
		height = size.height()
		width = size.width()

		#out_height = int(height * 0.4)
		#plot_height = height - out_height
		#docks = [self.output_dock, self.plot_dock]
		#sizes = [out_height, plot_height]
		#self.resizeDocks(docks, sizes, Qt.Vertical)

	def create_actions(self):
		self.open_act = QAction("Open Project", self)
		self.open_act.setShortcut(QKeySequence.Open)
		self.open_act.triggered.connect(self.on_open_project)

		self.save_act = QAction("Save Project", self)
		self.save_act.setShortcut(QKeySequence.Save)
		self.save_act.triggered.connect(self.on_save_project)

		self.saveas_act = QAction("Save Project As", self)
		self.saveas_act.setShortcut(QKeySequence.Open)
		self.saveas_act.triggered.connect(self.on_save_project_as)

		self.close_act = QAction("Close Project", self)
		self.close_act.setShortcut(QKeySequence.Close)
		self.close_act.triggered.connect(self.on_close_project)

		self.import_count_act = QAction("Read Counts Matrix", self)
		self.import_count_act.triggered.connect(self.on_import_read_counts)

		self.import_tpm_act = QAction("TPM Matrix")
		self.import_tpm_act.triggered.connect(self.on_import_tpm_matrix)

		self.import_fpkm_act = QAction("FPKM Matrix")
		self.import_fpkm_act.triggered.connect(self.on_import_fpkm_matrix)

		self.import_salmon_act = QAction("From Salmon")
		self.import_salmon_act.triggered.connect(self.on_import_from_salmon)

		self.import_sailfish_act = QAction("From Sailfish")
		self.import_sailfish_act.triggered.connect(self.on_import_from_sailfish)

		self.import_kallisto_act = QAction("From Kallisto")
		self.import_kallisto_act.triggered.connect(self.on_import_from_kallisto)

		self.import_rsem_act = QAction("From RSEM")
		self.import_rsem_act.triggered.connect(self.on_import_from_rsem)

		self.import_sample_act = QAction("Import Sample Information", self)
		self.import_sample_act.triggered.connect(self.on_import_sample_information)

		self.import_deg_act = QAction("Import DEG List", self)
		self.import_deg_act.triggered.connect(self.on_import_degs_list)

		self.ingname_act = QAction("Gene Annotation", self)
		self.ingname_act.triggered.connect(self.on_import_gene_annotation)

		self.ingoann_act = QAction("GO Annotation", self)
		self.ingoann_act.triggered.connect(self.on_import_go_annotation)

		self.inkegga_act = QAction("KEGG Annotation", self)
		self.inkegga_act.triggered.connect(self.on_import_kegg_annotation)

		self.input_dock_act = self.input_dock.toggleViewAction()
		self.input_dock_act.setText("Show Input Tables")
		self.output_dock_act = self.output_dock.toggleViewAction()
		self.output_dock_act.setText("Show Output Tables")
		#self.plot_dock_act = self.plot_dock.toggleViewAction()
		#self.plot_dock_act.setText("Show Output Plots")

		self.deseq_degs_act = QAction("Identify DEGs by DESeq2", self)
		self.deseq_degs_act.triggered.connect(self.do_identify_degs_by_deseq2)

		self.edger_degs_act = QAction("Identify DEGs by edgeR", self)
		self.edger_degs_act.triggered.connect(self.do_identify_degs_by_edger)

		self.show_degs_act = QAction("Extract Identified DEGs", self)
		self.show_degs_act.triggered.connect(self.do_extract_identifid_degs)

		self.sample_cluster_act = QAction("Sample Cluster Analysis", self)

		self.kegg_enrich_act = QAction("KEGG Enrichment")
		self.go_enrich_act = QAction("GO Enrichment")

		self.dist_plot_act = QAction("Distribution Plot", self)
		self.dist_plot_act.triggered.connect(self.do_plot_degs_dist)

		self.venn_plot_act = QAction("Venn Diagram", self)
		self.venn_plot_act.triggered.connect(self.do_plot_degs_venn)

		self.upset_plot_act = QAction("Upset Plot", self)
		self.upset_plot_act.triggered.connect(self.do_plot_degs_upset)

		self.volcano_plot_act = QAction("Volcano Plot", self)
		self.volcano_plot_act.triggered.connect(self.do_plot_degs_volcano)

		self.ma_plot_act = QAction("MA plot", self)
		self.ma_plot_act.triggered.connect(self.do_plot_degs_ma)

		self.exit_act = QAction("Exit", self)
		self.exit_act.setShortcut(QKeySequence.Quit)
		self.exit_act.triggered.connect(self.close)

		self.global_set_act = QAction("Global Settings")
		self.global_set_act.setShortcut(QKeySequence.Preferences)
		self.global_set_act.triggered.connect(self.do_open_setting_dialog)

		self.install_act = QAction("R Package Manager", self)
		self.install_act.triggered.connect(self.on_open_package_dialog)

		self.about_act = QAction("About", self)
		self.about_act.triggered.connect(self.on_open_about_dialog)

		#tool bar actions
		self.deg_act = QAction(QIcon("icons/deg.svg"), "Identify Differentially Expressed Genes", self)
		self.cluster_act = QAction(QIcon("icons/cluster.svg"), "Cluster Analysis", self)
		self.enrich_act = QAction(QIcon("icons/enrich.svg"), "Enrichment Analysis", self)

	def create_menus(self):
		self.file_menu = self.menuBar().addMenu("&File")
		self.file_menu.addAction(self.open_act)
		self.file_menu.addAction(self.save_act)
		self.file_menu.addAction(self.saveas_act)
		self.file_menu.addAction(self.close_act)
		self.file_menu.addSeparator()
		self.exp_menu = self.file_menu.addMenu("Import Expression Data")
		self.exp_menu.addAction(self.import_count_act)
		self.exp_menu.addAction(self.import_tpm_act)
		self.exp_menu.addAction(self.import_fpkm_act)
		self.exp_menu.addSeparator()
		self.exp_menu.addAction(self.import_salmon_act)
		self.exp_menu.addAction(self.import_sailfish_act)
		self.exp_menu.addAction(self.import_kallisto_act)
		self.exp_menu.addAction(self.import_rsem_act)

		self.file_menu.addAction(self.import_sample_act)
		self.file_menu.addAction(self.import_deg_act)

		self.import_menu = self.file_menu.addMenu("&Import Annotation Data")
		self.import_menu.addAction(self.ingname_act)
		self.import_menu.addAction(self.ingoann_act)
		self.import_menu.addAction(self.inkegga_act)
		self.file_menu.addSeparator()
		self.file_menu.addAction(self.exit_act)

		self.edit_menu = self.menuBar().addMenu("&Edit")
		#self.edit_menu.addAction(self.global_set_act)

		self.view_menu = self.menuBar().addMenu("&View")
		self.view_menu.addAction(self.input_dock_act)
		self.view_menu.addAction(self.output_dock_act)
		#self.view_menu.addAction(self.plot_dock_act)

		self.anal_menu = self.menuBar().addMenu("&Analysis")
		deg_menu = self.anal_menu.addMenu("DEGs Analysis")
		deg_menu.addAction(self.deseq_degs_act)
		deg_menu.addAction(self.edger_degs_act)
		deg_menu.addSeparator()
		deg_menu.addAction(self.show_degs_act)
		self.anal_menu.addAction(self.sample_cluster_act)
		self.enrich_menu = self.anal_menu.addMenu("Enrichment Analysis")
		self.enrich_menu.addAction(self.go_enrich_act)
		self.enrich_menu.addAction(self.kegg_enrich_act)

		self.plot_menu = self.menuBar().addMenu("&Plots")
		deg_plots = self.plot_menu.addMenu("DEGs")
		deg_plots.addAction(self.dist_plot_act)
		deg_plots.addAction(self.venn_plot_act)
		deg_plots.addAction(self.upset_plot_act)
		deg_plots.addSeparator()
		deg_plots.addAction(self.volcano_plot_act)
		deg_plots.addAction(self.ma_plot_act)
		
		self.tool_menu = self.menuBar().addMenu("&Tools")
		self.tool_menu.addAction(self.install_act)
		self.tool_menu.addSeparator()
		self.tool_menu.addAction(self.global_set_act)

		self.help_menu = self.menuBar().addMenu("&Help")
		self.help_menu.addAction(self.about_act)

	def create_toolbar(self):
		self.tool_bar = self.addToolBar("Show Toolbar")
		self.tool_bar.setIconSize(QSize(28, 28))
		self.tool_bar.setMovable(False)

		#self.tool_bar.addWidget(self.input)
		#self.tool_bar.addWidget(self.button)

		self.tool_bar.addAction(self.deg_act)
		self.tool_bar.addAction(self.cluster_act)
		self.tool_bar.addAction(self.enrich_act)

		spacer = RNASuiteSpacerWidget(self)
		self.tool_bar.addWidget(spacer)

		self.wait_spinner = RNASuiteWaitingSpinner(self)
		self.run_progress = self.tool_bar.addWidget(self.wait_spinner)
		self.run_progress.setVisible(False)

	def create_statusbar(self):
		self.status_bar = self.statusBar()

	@Slot()
	def on_open_package_dialog(self):
		dlg = RNASuitePackageManagerDialog(self)
		dlg.exec()

	@Slot()
	def do_open_setting_dialog(self):
		dlg = RNASuiteGlobalSettingDialog(self)
		dlg.exec()

	def import_data_file(self, title, reader):
		delimiter = None
		file, _ = QFileDialog.getOpenFileName(self,
			caption = title,
			filter = "Table files (*.csv *.tsv *.xls *.xlsx *.txt);;All files (*.*)"
		)

		if not file:
			return

		if not file.endswith(('.csv', '.tsv', '.xls', '.xlsx')):
			delimiter = RNASuiteColumnSeparatorDialog.get_delimiter(self)

		match reader:
			case 'read_count':
				self.input_list.import_read_count(file, delimiter)

			case 'sample_info':
				self.input_list.import_sample_info(file, delimiter)





	@Slot()
	def on_import_read_counts(self):
		self.import_data_file("Select Read Count Matrix File", 'read_count')

	@Slot()
	def on_import_fpkm_matrix(self):
		self.import_data_file("Select FPKM Matrix File", 'gene_fpkm')

	@Slot()
	def on_import_tpm_matrix(self):
		self.import_data_file("Select TPM Matrix File", 'gene_tpm')

	@Slot()
	def on_import_degs_list(self):
		self.import_data_file("Select DEG List File", 'degs_list')

	@Slot()
	def on_import_sample_information(self):
		self.import_data_file("Select sample information file", 'sample_info')

	@Slot()
	def on_import_gene_annotation(self):
		self.import_data_file("Select gene annotation file", 'gene_annot')

	@Slot()
	def on_import_go_annotation(self):
		self.import_data_file("Select GO annotation file", 'go_annot')

	@Slot()
	def on_import_kegg_annotation(self):
		self.import_data_file("Select KEGG annotation file", 'kegg_annot')

	@Slot()
	def on_import_from_dir(self, quant_name=None, count_col=None, fpkm_col=None, tpm_col=None):
		folder = QFileDialog.getExistingDirectory(self,
			caption = "Select a Directory",
			options = QFileDialog.ShowDirOnly | QFileDialog.DontResolveSymlinks
		)

		if not folder:
			return

		quant_files = QDirIterator(folder, QDir.Files, QDirIterator.Subdirectories)

		while quant_files.hasNext():
			quant_file = quant_files.next()
			file_info = QFileInfo(quant_file)
			file_name = file_info.fileName()

			if file_name == quant_name:
				sample_id = QFileInfo(file_info.path()).baseName()

	@Slot()
	def on_import_from_rsem(self):
		pass

	@Slot()
	def on_import_from_salmon(self):
		pass

	@Slot()
	def on_import_from_kallisto(self):
		pass

	@Slot()
	def on_import_from_sailfish(self):
		pass

	def run_analysis_worker(self, worker):
		worker.error.connect(self.show_error_message)
		worker.start()

	def has_running_worker(self):
		if self.wait_spinner.running:
			QMessageBox.warning(self, "Warning", "A task is already running")
			return True

		return False

	def has_none_deg_data(self):
		if self.table_widgets.is_empty('read_count'):
			self.main_error.prompt('count_error')
			return True

		if self.table_widgets.is_empty('sample_info'):
			self.main_error.prompt('sample_error')
			return True

		return False

	@Slot()
	def do_identify_degs_by_deseq2(self):
		if self.has_running_worker():
			return

		#if self.has_none_deg_data():
		#	return

		#defines = self.stored_params.get('degs', {})
		#samples = self.table_widgets.get_data('sample_info')
		#samples = self.input_list.sample_info
		#dataset = {c: list(samples[c].unique()) for c in samples.columns}

		params = RNASuiteDeseqParameterDialog.get_params(self)

		if not params:
			return

		#read_counts = self.input_list.get_tight('read_counts')
		#sample_info = self.input_list.get_tight('sample_info')
		#params['counts'] = self.input_list.read_counts
		#params['samples'] = self.input_list.sample_info

		worker = RNASuiteDeseqIdentifyWorker(self, params)
		self.run_analysis_worker(worker)
		#params['tool'] = 'deseq'
		#self.stored_params['degs'] = params

	@Slot()
	def do_identify_degs_by_edger(self):
		if self.has_running_worker():
			return

		if self.has_none_deg_data():
			return

		defines = self.stored_params.get('degs', {})
		samples = self.table_widgets.get_data('sample_info')
		dataset = {c: list(samples[c].unique()) for c in samples.columns}
		params = RNASuiteEdgerParameterDialog.get_params(self, defines, dataset)

		if not params:
			return

		read_counts = self.table_widgets.get_tight('read_count')
		sample_info = self.table_widgets.get_tight('sample_info')
		worker = RNASuiteEdgerDEGWorker(self, read_counts, sample_info, params)
		self.run_analysis_worker(worker)
		params['tool'] = 'edger'
		self.stored_params['degs'] = params

	def has_none_identified_degs(self):
		if 'degs' not in self.stored_params or not self.table_widgets.has_table('degs_list'):
			self.main_error.prompt('degs_error')
			return True

		return False

	@Slot()
	def do_extract_identifid_degs(self):
		params = RNASuiteExtractDegsParameterDialog.get_params(self)

		if not params:
			return

		worker = RNASuiteDegsExtractWorker(self, params)
		self.run_analysis_worker(worker)

	@Slot()
	def do_plot_degs_dist(self):
		params = RNASuiteDegsDistPlotParameterDialog.get_params(self)

		if not params:
			return

		worker = RNASuiteDegsDistPlotWorker(self, params)
		self.run_analysis_worker(worker)


	@Slot()
	def do_plot_degs_venn(self):
		params = RNASuiteDegsVennPlotParameterDialog.get_params(self)

		if not params:
			return

		worker = RNASuiteDegsVennPlotWorker(self, params)
		self.run_analysis_worker(worker)

	@Slot()
	def do_plot_degs_upset(self):
		params = RNASuiteDegsUpsetPlotParameterDialog.get_params(self)

		if not params:
			return

		worker = RNASuiteDegsUpsetPlotWorker(self, params)
		self.run_analysis_worker(worker)

	@Slot()
	def do_plot_degs_volcano(self):
		params = RNASuiteDegsVolcanoPlotParameterDialog.get_params(self)
		
		if not params:
			return

		worker = RNASuiteDegsVolcanoPlotWorker(self, params)
		self.run_analysis_worker(worker)

	@Slot()
	def do_plot_degs_ma(self):
		pass

	@Slot()
	def on_open_about_dialog(self):
		about_str = "an interactive tool for identifying differentailly expressed genes and enrichment analysis"
		QMessageBox.about(self, "About", about_str)

	def create_input_dock(self):
		#self.input_tabs = QTabWidget(self)
		self.input_list = RNASuiteInputListWidget(self)
		self.input_list.show_table.connect(self._on_show_table_data)
		self.input_dock = QDockWidget("Input", self)
		#self.input_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
		self.input_dock.setAllowedAreas(Qt.LeftDockWidgetArea)
		self.input_dock.setWidget(self.input_list)
		self.addDockWidget(Qt.LeftDockWidgetArea, self.input_dock)
		#self.setCentralWidget(self.input_tabs)

	def create_output_dock(self):
		#self.output_tabs = QTabWidget(self)
		self.output_list = RNASuiteOutputTreeWidget(self)
		self.output_list.show_table.connect(self._on_show_table_data)
		self.output_list.show_plot.connect(self._on_show_table_plot)
		self.output_list.remove_plot.connect(self.plot_viewer.remove_plot)
		self.output_list.show_panel.connect(self.plot_viewer.show_panel)
		self.output_dock = QDockWidget("Output", self)
		#self.output_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
		self.output_dock.setAllowedAreas(Qt.LeftDockWidgetArea)
		self.output_dock.setWidget(self.output_list)
		self.addDockWidget(Qt.LeftDockWidgetArea, self.output_dock)

	def create_main_widgets(self):
		self.data_table = RNASuitePandasTable(self)
		self.plot_viewer = RNASuitePlotViewer(self)
		self.plot_viewer.error.connect(self.show_error_message)
		#self.main_tabs = QTabWidget(self)
		#self.main_tabs.addTab(self.data_table, "Table")
		#self.main_tabs.addTab(self.plot_viewer, "Plot")
		self.stack_widget = QStackedWidget(self)
		self.stack_widget.addWidget(self.data_table)
		self.stack_widget.addWidget(self.plot_viewer)
		self.stack_widget.setCurrentIndex(1)
		self.setCentralWidget(self.stack_widget)
		#self.setCentralWidget(self.main_tabs)

	#def create_plot_dock(self):
		#self.plot_stack = RNASuitePlotStackedWidget(self)
		#self.plot_dock = QDockWidget("Plot", self)
		#self.plot_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
		#self.plot_dock.setAllowedAreas(Qt.RightDockWidgetArea)
		#self.plot_dock.setWidget(self.plot_stack)
		#self.addDockWidget(Qt.RightDockWidgetArea, self.plot_dock)
		

	@Slot()
	def _on_show_table_data(self, data):
		#dlg = RNASuiteShowPandasDataDialog(self, data)
		#dlg.exec()
		self.data_table.update_data(data)
		self.stack_widget.setCurrentIndex(0)

	@Slot()
	def _on_show_table_plot(self, plot_id):
		self.plot_viewer.show_plot(plot_id)
		self.stack_widget.setCurrentIndex(1)

	@Slot(str)
	def show_warn_message(self, warn):
		QMessageBox.warning(self, "Warning", warn)

	@Slot(str)
	def show_error_message(self, error):
		QMessageBox.critical(self, "Error occurred", error)

	@Slot(str)
	def show_status_message(self, text):
		self.status_bar.showMessage(text)

