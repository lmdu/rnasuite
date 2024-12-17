import os
import math

import pandas
from rchitect.utils import Rhome

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *
from PySide6.QtSvgWidgets import *

from utils import *
from models import *
from config import *
from threads import *

__all__ = [
	'create_parameter_widget',
	'set_parameter_widget_value',
	'get_parameter_widget_value',
	'RNASuiteInputListWidget',
	'RNASuiteOutputTreeWidget',
	'RNASuitePackageInstallButton',
	'RNASuiteWaitingSpinner',
	'RNASuitePackageTreeView',
	'RNASuitePackageInstallMessage',
	'RNASuiteSpacerWidget',
	'RNASuiteMultipleSelect',
	'RNASuiteRGeneralSettingPage',
	'RNASuiteColorButton',
	'RNASuiteContrastVersusWidget',
	'RNASuiteColorGroups',
	'RNASuiteAccordionWidget',
	'RNASuiteAxisLimitWidget'
]

def create_parameter_widget(param):
	match param.type:
		case 'int':
			widget = QSpinBox()
			widget.setRange(*param.range)
			widget.setSingleStep(param.step)

		case 'float':
			widget = QDoubleSpinBox()
			widget.setRange(*param.range)
			widget.setSingleStep(param.step)
			widget.setDecimals(param.decimal)

		case 'str':
			widget = QLineEdit()

		case 'list':
			widget = QComboBox()
			widget.addItems(param.options)

		case 'bool':
			widget = QCheckBox()

		case 'select':
			widget = RNASuiteMultipleSelect()
			widget.add_items(param.options)

		case 'text':
			widget = QPlainTextEdit()

		case 'color':
			widget = RNASuiteColorButton()

		case 'colors':
			widget = RNASuiteColorGroups()

		case 'contrast':
			widget = RNASuiteContrastVersusWidget()

		case 'limit':
			widget = RNASuiteAxisLimitWidget()
			widget.set_ranges(*param.range)
			widget.set_steps(param.step)
			widget.set_decimals(param.decimal)

	return widget

def set_parameter_widget_value(widget, value=None, index=False):
	if value is None:
		return

	match widget:
		case QAbstractSpinBox():
			widget.setValue(value)

		case QLineEdit():
			widget.setText(value)

		case QComboBox():
			if index:
				widget.setCurrentIndex(value)
			else:
				widget.setCurrentText(value)

		case QCheckBox():
			widget.setChecked(value)

		case QPlainTextEdit():
			widget.setPlainText(value)

		case RNASuiteMultipleSelect():
			widget.set_text(value)

		case RNASuiteColorButton():
			widget.set_color(value)

		case RNASuiteColorGroups():
			widget.set_colors(value)

		case RNASuiteContrastVersusWidget():
			widget.set_contrasts(value)

		case RNASuiteAxisLimitWidget():
			widget.set_limits(value)

def get_parameter_widget_value(widget, index=False):
	match widget:
		case QAbstractSpinBox():
			value = widget.value()

		case QLineEdit():
			value = widget.text()

		case QComboBox():
			if index:
				value = widget.currentIndex()
			else:
				value = widget.currentText()

		case QCheckBox():
			value = widget.isChecked()

		case QPlainTextEdit():
			value = widget.toPlainText()

		case RNASuiteMultipleSelect():
			value = widget.get_text()

		case RNASuiteColorButton():
			value = widget.get_color()

		case RNASuiteColorGroups():
			value = widget.get_colors()

		case RNASuiteContrastVersusWidget():
			value = widget.get_contrasts()

		case RNASuiteAxisLimitWidget():
			value = widget.get_limits()

	return value

class RNASuiteAxisLimitWidget(QWidget):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.create_widgets()
		self.set_layouts()

	def sizeHint(self):
		return QSize(100, 10)

	def create_widgets(self):
		self.min_limit = QDoubleSpinBox(self)
		self.min_limit.setDecimals(2)
		self.join_label = QLabel("~", self)
		self.max_limit = QDoubleSpinBox(self)
		self.max_limit.setDecimals(2)

	def set_layouts(self):
		main_layout = QHBoxLayout()
		main_layout.setContentsMargins(0, 0, 0, 0)
		main_layout.addWidget(self.min_limit, 1)
		main_layout.addWidget(self.join_label)
		main_layout.addWidget(self.max_limit, 1)
		self.setLayout(main_layout)

	def set_ranges(self, minv, maxv):
		self.min_limit.setRange(minv, maxv)
		self.max_limit.setRange(minv, maxv)

	def set_steps(self, step=1):
		self.min_limit.setSingleStep(step)
		self.max_limit.setSingleStep(step)

	def set_decimals(self, decimal=1):
		self.min_limit.setDecimals(decimal)
		self.max_limit.setDecimals(decimal)

	def get_limits(self):
		return [self.min_limit.value(), self.max_limit.value()]

	def set_limits(self, limits):
		self.min_limit.setValue(limits[0])
		self.max_limit.setValue(limits[1])

class RNASuiteAccordionHeader(QPushButton):
	def __init__(self, parent=None, title=None):
		super().__init__(parent)
		self.setText(title)
		self.setCheckable(True)
		self.expand_icon = QIcon('icons/down.svg')
		self.collapse_icon = QIcon('icons/right.svg')
		self.setIcon(self.collapse_icon)
		self.clicked.connect(self._on_clicked)

	def _on_clicked(self):
		if self.isChecked():
			self.setIcon(self.expand_icon)

		else:
			self.setIcon(self.collapse_icon)

class RNASuiteAccordionItem(QWidget):
	def __init__(self, parent=None, title=None):
		super().__init__(parent)
		
		self.header = RNASuiteAccordionHeader(self, title)
		self.content = QWidget(self)
		self.content.setVisible(False)
		self.header.toggled.connect(self.content.setVisible)

		self.set_layouts()

	def sizeHint(self):
		return QSize(100, 10)

	def set_layouts(self):
		self.main_layout = QVBoxLayout()
		self.main_layout.setContentsMargins(0, 0, 0, 0)
		self.main_layout.addWidget(self.header)
		self.main_layout.addWidget(self.content)
		self.setLayout(self.main_layout)

		self.content_layout = QVBoxLayout()
		self.content_layout.setContentsMargins(10, 0, 0, 10)
		self.content.setLayout(self.content_layout)

	def add_widgets(self, widgets):
		for label, widget in widgets:
			if label:
				self.content_layout.addWidget(label)

			self.content_layout.addWidget(widget)

class RNASuiteAccordionWidget(QScrollArea):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.setFrameStyle(QFrame.NoFrame)
		#self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
		self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		self.setWidgetResizable(True)

		self.main_widget = QWidget(self)
		self.main_layout = QVBoxLayout()
		self.main_layout.setContentsMargins(0, 0, 2, 0)
		self.main_widget.setLayout(self.main_layout)
		self.setWidget(self.main_widget)

	def sizeHint(self):
		return QSize(100, 10)

	def add_accordions(self, title, widgets):
		accordion = RNASuiteAccordionItem(self, title)
		accordion.add_widgets(widgets)
		self.main_layout.addWidget(accordion)

	def add_stretcher(self):
		self.main_layout.addStretch()

class RNASuiteInputListItem(QWidget):
	def __init__(self, parent=None, title=None, content=None, meta=None):
		super().__init__(parent)

		title_label = QLabel("<b><small>{}</small></b>".format(title), self)
		content_label = QLabel("<small>{}</small>".format(content), self)
		meta_label = QLabel("<small>{}</small>".format(meta))

		main_layout = QVBoxLayout()
		main_layout.setSpacing(0)
		main_layout.setContentsMargins(0, 0, 0, 0)
		main_layout.addWidget(title_label)
		main_layout.addWidget(content_label)
		main_layout.addWidget(meta_label, alignment=Qt.AlignRight)
		self.setLayout(main_layout)

class RNASuiteInputListWidget(QListWidget):
	show_table = Signal(pandas.DataFrame)

	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent
		self.setIconSize(QSize(28, 28))
		self.setSpacing(2)
		self.itemClicked.connect(self._on_open_table)

	def sizeHint(self):
		return QSize(200, 200)

	def add_item_widget(self, index, icon, title, text, meta):
		widget = RNASuiteInputListItem(self, title, text, meta)
		item = QListWidgetItem()
		item.setIcon(QIcon(icon))
		self.insertItem(index, item)
		self.setItemWidget(item, widget)
		item.setSizeHint(widget.sizeHint())
		return item

	def read_file(self, file, delimiter=None):
		if file.endswith('.csv'):
			data = pandas.read_csv(file, index_col=0)

		elif file.endswith('.tsv'):
			data = pandas.read_csv(file, index_col=0, sep='\t')

		elif file.endswith(('.xls', '.xlsx')):
			data = pandas.read_excel(file, index_col=0)

		else:
			data = pandas.read_table(file, index_col=0, sep=delimiter)

		return data

	def has_data(self, dtype):
		if hasattr(self, dtype):
			x = getattr(self, dtype)

			if len(x) > 0:
				return True

		return False

	def get_groups(self):
		return {c: list(self.sample_info[c].unique()) for c in self.sample_info.columns}

	def get_tight(self, dtype):
		match dtype:
			case 'read_counts':
				return self.read_counts.to_dict(orient='tight')

			case 'sample_info':
				return self.sample_info.to_dict(orient='tight')

	def import_read_count(self, file, delimiter):
		self.read_counts = self.read_file(file, delimiter)
		self.count_file = file
		genes = len(self.read_counts)
		samples = len(self.read_counts.columns)

		self.count_item = self.add_item_widget(0,
			'icons/c.svg',
			'Read Counts',
			os.path.basename(file),
			'Genes: {} Samples: {}'.format(genes, samples)
		)

	def import_sample_info(self, file, delimiter):
		self.sample_info = self.read_file(file, delimiter)
		self.sample_file = file
		samples = len(self.sample_info)

		self.sample_item = self.add_item_widget(1,
			'icons/s.svg',
			'Sample Information',
			os.path.basename(file),
			'Samples: {}'.format(samples)
		)

	@Slot()
	def _on_open_table(self, item):
		match item:
			case self.count_item:
				self.show_table.emit(self.read_counts)

			case self.sample_item:
				self.show_table.emit(self.sample_info)

class RNASuiteOutputTreeWidget(QTreeView):
	show_table = Signal(object)
	show_panel = Signal(str, int, int)
	show_plot = Signal(int)
	remove_plot = Signal(int)

	def __init__(self, parent=None):
		super().__init__(parent)

		self.setRootIsDecorated(False)
		self.clicked.connect(self._on_row_clicked)
		self.create_model()

	def sizeHint(self):
		return QSize(200, 500)

	def create_model(self):
		self._model = RNASuiteOutputTreeModel(self)
		self.setModel(self._model)
		self.header().setStretchLastSection(False)
		self.header().setSectionResizeMode(0, QHeaderView.Stretch)
		self.header().setSectionResizeMode(1, QHeaderView.ResizeToContents)

	@Slot()
	def on_receive_data(self, results):
		if not results:
			return

		for result in results:
			output = self._model.get_output_by_name(result['name'])

			if output:
				self._model.update_output(result, output.id)

				if output.type == 'plot':
					self.remove_plot.emit(output.code)
			else:
				self._model.add_output(result)

	def _on_row_clicked(self, index):
		row = index.row()
		plot = self._data.iloc[row, 2]
		did = self._data.iloc[row, 3]
		ptype = self._data.iloc[row, 4]
		pyid = self._data.iloc[row, 5]

		if plot:
			self.show_plot.emit(did)
			self.show_panel.emit(ptype, did, pyid)
			self._data.loc[self._data['plot'] == 1, 'update'] = 0
			self._data.iloc[row, 1] = 1

		else:
			data = self.datasets[did]
			self.show_table.emit(data)
			self._data.iloc[row, 1] = 0

		self._model.load_data(self._data)

class RNASuiteSpacerWidget(QWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

#from https://www.pythonguis.com/widgets/qcolorbutton-a-color-selector-tool-for-pyqt/
class RNASuiteColorButton(QPushButton):
	color_changed = Signal(tuple)

	def __init__(self, parent=None):
		super().__init__(parent)

		self._color = None
		self._default = '#ffffff'
		self.pressed.connect(self.on_select_color)
		self.set_color(self._default)
		self.setFixedSize(QSize(28, 16))

	def set_color(self, color):
		if color != self._color:
			self._color = color
			self.color_changed.emit(QColor(color).getRgbF()[0:3])

		if self._color:
			self.setStyleSheet("background-color: {}".format(self._color))

		else:
			self.setStyleSheet("")

	def get_color(self):
		return self._color

	@Slot()
	def on_select_color(self):
		dlg = QColorDialog()

		if self._color:
			dlg.setCurrentColor(QColor(self._color))

		if dlg.exec():
			self.set_color(dlg.currentColor().name())

	def mousePressEvent(self, event):
		if event.button() == Qt.RightButton:
			self.set_color(self._default)

		return super().mousePressEvent(event)

#from https://github.com/z3ntu/QtWaitingSpinner
class RNASuiteWaitingSpinner(QWidget):
	def __init__(self, parent):
		super().__init__(parent)

		# WAS IN initialize()
		self._color = QColor(Qt.black)
		self._roundness = 70.0
		self._minimumTrailOpacity = 3.14159265358979323846
		self._trailFadePercentage = 70.0
		self._revolutionsPerSecond = 1.57079632679489661923
		self._numberOfLines = 10
		self._lineLength = 6
		self._lineWidth = 3
		self._innerRadius = 4
		self._currentCounter = 0
		self._isSpinning = False

		self._timer = QTimer(self)
		self._timer.timeout.connect(self.rotate)
		self.updateSize()
		self.updateTimer()
		#self.setVisible(False)
		# END initialize()

		self.setAttribute(Qt.WA_TranslucentBackground)

	def paintEvent(self, QPaintEvent):
		painter = QPainter(self)
		painter.fillRect(self.rect(), Qt.transparent)
		painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

		if self._currentCounter >= self._numberOfLines:
			self._currentCounter = 0

		painter.setPen(Qt.NoPen)
		for i in range(0, self._numberOfLines):
			painter.save()
			painter.translate(self._innerRadius + self._lineLength, self._innerRadius + self._lineLength)
			rotateAngle = float(360 * i) / float(self._numberOfLines)
			painter.rotate(rotateAngle)
			painter.translate(self._innerRadius, 0)
			distance = self.lineCountDistanceFromPrimary(i, self._currentCounter, self._numberOfLines)
			color = self.currentLineColor(distance, self._numberOfLines, self._trailFadePercentage,
										  self._minimumTrailOpacity, self._color)
			painter.setBrush(color)
			rect = QRect(0, int(-self._lineWidth / 2), int(self._lineLength), int(self._lineWidth))
			painter.drawRoundedRect(rect, self._roundness, self._roundness, Qt.SizeMode.RelativeSize)
			painter.restore()

	@property
	def running(self):
		return self._isSpinning

	def toggle(self, state):
		if state:
			self.start()
		else:
			self.stop()

	def start(self):
		self._isSpinning = True
		#self.setVisible(True)

		if not self._timer.isActive():
			self._timer.start()
			self._currentCounter = 0

	def stop(self):
		self._isSpinning = False
		#self.setVisible(False)

		if self._timer.isActive():
			self._timer.stop()
			self._currentCounter = 0

	def rotate(self):
		self._currentCounter += 1
		if self._currentCounter >= self._numberOfLines:
			self._currentCounter = 0
		self.update()

	def updateSize(self):
		size = int((self._innerRadius + self._lineLength) * 2)
		self.setFixedSize(size, size)

	def updateTimer(self):
		self._timer.setInterval(int(1000 / (self._numberOfLines * self._revolutionsPerSecond)))

	def lineCountDistanceFromPrimary(self, current, primary, totalNrOfLines):
		distance = primary - current
		if distance < 0:
			distance += totalNrOfLines
		return distance

	def currentLineColor(self, countDistance, totalNrOfLines, trailFadePerc, minOpacity, colorinput):
		color = QColor(colorinput)
		if countDistance == 0:
			return color
		minAlphaF = minOpacity / 100.0
		distanceThreshold = int(math.ceil((totalNrOfLines - 1) * trailFadePerc / 100.0))
		if countDistance > distanceThreshold:
			color.setAlphaF(minAlphaF)
		else:
			alphaDiff = color.alphaF() - minAlphaF
			gradient = alphaDiff / float(distanceThreshold + 1)
			resultAlpha = color.alphaF() - gradient * countDistance
			# If alpha is out of bounds, clip it.
			resultAlpha = min(1.0, max(0.0, resultAlpha))
			color.setAlphaF(resultAlpha)
		return color

class RNASuitePackageInstallButton(QPushButton):
	install = Signal(QWidget, str, str)

	def __init__(self, package, repository, parent=None, action='Install'):
		super().__init__(parent)
		self.setText(action)
		self.clicked.connect(self.on_clicked)

		self.action = action
		self.package = package
		self.repository = repository

	@Slot()
	def on_clicked(self):
		self.install.emit(self, self.package, self.repository)

class RNASuitePackageTreeView(QTreeView):
	error = Signal(str)
	message = Signal(str)
	started = Signal()
	stopped = Signal()

	def __init__(self, parent=None):
		super().__init__(parent)
		self.header().setStretchLastSection(False)

		self.get_rscript()
		self.create_model()
		self.update_tree()

		self.task = None

	def sizeHint(self):
		return QSize(750, 400)

	def on_error_occurred(self, msg):
		self.error.emit(str(msg))

	def get_rscript(self):
		if os.name == 'nt':
			self.rscript = os.path.join(Rhome(), 'bin', 'Rscript.exe')
		else:
			self.rscript = os.path.join(Rhome(), 'bin', 'Rscript')

		if not os.path.isfile(self.rscript):
			self.on_error_occurred("Could not find Rscript: {}".format(self.rscript))

	def update_colum_width(self):
		self.header().setSectionResizeMode(0, QHeaderView.ResizeToContents)
		self.header().setSectionResizeMode(1, QHeaderView.Stretch)
		self.header().setSectionResizeMode(2, QHeaderView.Interactive)
		self.header().setSectionResizeMode(3, QHeaderView.Interactive)
		self.header().setSectionResizeMode(4, QHeaderView.Interactive)
		self.expandAll()
		self.setWordWrap(True)

	def create_model(self):
		self._model = QStandardItemModel(self)
		self.setModel(self._model)

	def update_tree(self):
		thread = RNASuitePackageVersionThread(self)
		thread.started.connect(self.on_task_started)
		thread.finished.connect(self.on_update_finished)
		thread.result.connect(self.update_data)
		thread.start()

	def update_data(self, versions):
		self._model.setHorizontalHeaderLabels([
			'Name', 'Description', 'Required', 'Version', 'Status'
		])

		self.packages = RNASuitePackages()

		for c in self.packages.get_order():
			for p in self.packages[c]:
				p.version = versions.get(p.name.lower(), '')
				p.status = compare_version(p.required, p.version)

		root = self._model.invisibleRootItem()
		for c in self.packages.get_order():
			parent = QStandardItem(c)
			font = QFont()
			font.setBold(True)
			parent.setFont(font)
			root.appendRow(parent)

			for p in self.packages[c]:
				if p.status == 0:
					first = QStandardItem(QIcon('icons/no.svg'), p.name)
					last = QStandardItem('install')
					btn = RNASuitePackageInstallButton(p.name, p.repository, self)
					btn.install.connect(self.on_install_package)
				elif p.status == -1:
					first = QStandardItem(QIcon('icons/up.svg'), p.name)
					last = QStandardItem('update')
					btn = RNASuitePackageInstallButton(p.name, p.repository, self, 'Update')
					btn.install.connect(self.on_install_package)
				else:
					first = QStandardItem(QIcon('icons/ok.svg'), p.name)
					last = QStandardItem('OK')
					btn = None

				row = [
					first,
					QStandardItem(p.description),
					QStandardItem(p.required),
					QStandardItem(p.version),
					last
				]
				parent.appendRow(row)

				if btn is not None:
					self.setIndexWidget(last.index(), btn)

		self.update_colum_width()

	@Slot()
	def on_task_started(self):
		self.started.emit()

	@Slot()
	def on_update_finished(self):
		self.stopped.emit()

	@Slot()
	def on_task_finished(self, exit_code):
		self.stopped.emit()
		
		if exit_code == 0:
			self.message.emit("Installation completed!")
		else:
			self.message.emit("Installation failed with code: {}".format(exit_code))

	@Slot()
	def on_task_output_ready(self):
		out = self.task.readAllStandardOutput()
		self.message.emit(out.data().decode())

	@Slot()
	def on_task_error_ready(self):
		err = self.task.readAllStandardError()
		self.message.emit(err.data().decode())

	@Slot()
	def on_install_package(self, btn, package, repository):
		if self.task_running():
			self.error.emit("A pakcage installation is running")
			return

		settings = QSettings()
		default, convert = RNASUITE_SETTINGS['R']['cran_mirror']
		cran_mirror = settings.value('R/cran_mirror', default, convert)
		default, convert = RNASUITE_SETTINGS['R']['cran_mirror']
		bioc_mirror = settings.value('R/bioc_mirror', default, convert)

		rcode = None
		if package == 'R':
			QDesktopServices.openUrl("https://www.r-project.org/")
		elif repository == 'CRAN':
			rcode = (
				"options(repos=c(CRAN='{}'));"
				"install.packages('{}')"
			).format(cran_mirror, package)
		elif repository == 'Bioconductor':
			rcode = (
				"options(BioC_mirror='{}');"
				"if (!require('BiocManager',quietly=TRUE)){{install.packages('BiocManager')}};"
				"BiocManager::install('{}');"
			).format(bioc_mirror, package)

		if rcode is None:
			return

		btn.setDisabled(True)

		self.task = QProcess(self)
		#self.task.errorOccurred.connect(self.on_error_occurred)
		self.task.started.connect(self.on_task_started)
		self.task.finished.connect(self.on_task_finished)
		self.task.readyReadStandardOutput.connect(self.on_task_output_ready)
		self.task.readyReadStandardError.connect(self.on_task_error_ready)
		self.task.start(self.rscript, ['-e', rcode])		

	def task_running(self):
		if self.task:
			if self.task.state() == QProcess.Running:
				return True

		return False

	def stop_task(self):
		if self.task:
			if self.task.state() == QProcess.Running:
				self.task.kill()

	@Slot()
	def update_version(self):
		self._model.clear()
		self.update_tree()

class RNASuitePackageInstallMessage(QTextBrowser):
	def __init__(self, parent=None):
		super().__init__(parent)

	def sizeHint(self):
		return QSize(100, 25)

class RNASuiteMultipleSelect(QComboBox):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.setPlaceholderText('None')
		self.setEditable(True)
		self.lineEdit().setReadOnly(True)
		self._model = QStandardItemModel(self)
		self.setModel(self._model)
		self._model.itemChanged.connect(self.on_item_changed)

	def add_items(self, texts):
		self.clear()
		self.lineEdit().clear()

		for text in texts:
			item = QStandardItem(text)
			item.setCheckable(True)
			item.setSelectable(False)
			self._model.appendRow(item)

	@Slot()
	def on_item_changed(self, item):
		texts = self.lineEdit().text().split(' + ')

		if not texts[0]:
			texts = []

		if item.checkState() == Qt.Checked:
			if item.text() not in texts:
				texts.append(item.text())
		else:
			if item.text() in texts:
				texts.remove(item.text())

		self.lineEdit().setText(' + '.join(texts))

	def get_text(self):
		return self.lineEdit().text()

	def set_text(self, text):
		if not text:
			return

		self.lineEdit().setText(text)

		for val in text.split(';'):
			for item in self._model.findItems(val):
				item.setCheckState(Qt.Checked)

class RNASuiteBrowseLineEdit(QWidget):
	text_changed = Signal(str)

	def __init__(self, parent=None):
		super().__init__(parent)
	
		self.input = QLineEdit(self)
		self.input.setReadOnly(True)
		self.browse = QPushButton(self)
		self.browse.setFlat(True)
		self.browse.setIcon(QIcon("icons/folder.svg"))

		self.input.textChanged.connect(self.on_text_chanaged)
		self.browse.clicked.connect(self.select_file)

		layout = QHBoxLayout()
		layout.setSpacing(0)
		layout.setContentsMargins(0,0,0,0)
		layout.addWidget(self.input)
		layout.addWidget(self.browse)

		self.setLayout(layout)

	@Slot()
	def on_text_chanaged(self, text):
		self.text_changed.emit(text)

	@Slot()
	def select_file(self):
		pfile, _ = QFileDialog.getOpenFileName(self)

		if pfile:
			self.input.setText(pfile)

	def get_text(self):
		return self.input.text()

	def set_text(self, text):
		self.input.setText(text)

class RNASuiteGlobalSettingPage(QWidget):
	section = None

	def __init__(self, parent=None):
		super().__init__(parent)
		self.settings = QSettings()
		self.layout = QVBoxLayout()
		self.setLayout(self.layout)
		self.widgets = AttrDict()
		self.register_widgets()
		self.register_events()
		self.read_settings()

	@property
	def params(self):
		pass

	def register_widget(self, key, label, widget):
		self.widgets[key] = widget

		if label is not None:
			self.layout.addWidget(label)

		self.layout.addWidget(widget)

	def register_widgets(self):
		for p in self.params:
			if p.label is not None:
				label = QLabel(p.label, self)
			else:
				label = None

			if p.type == QFile:
				widget = RNASuiteBrowseLineEdit(self)
			elif p.type == str:
				widget = QLineEdit(self)
			elif p.type == list:
				widget = QListWidget(self)
			elif p.type == bool:
				widget = QCheckBox(p.text, self)
				
			self.register_widget(p.option, label, widget)

	def register_events(self):
		pass

	def read_settings(self):
		params = RNASUITE_SETTINGS[self.section]

		for param in params:
			default, convert = params[param]
			key = "{}/{}".format(self.section, param)
			value = self.settings.value(key, default, convert)
			widget = self.widgets[param]

			if isinstance(widget, QLineEdit):
				widget.setText(value)
			elif isinstance(widget, RNASuiteBrowseLineEdit):
				widget.set_text(value)

	def write_settings(self):
		params = RNASUITE_SETTINGS[self.section]

		for param in params:
			key = "{}/{}".format(self.section, param)
			widget = self.widgets[param]

			if isinstance(widget, QLineEdit):
				value = widget.text()
			elif isinstance(widget, RNASuiteBrowseLineEdit):
				value = widget.get_text()

			self.settings.setValue(key, value)

	def restore_settings(self):
		params = RNASUITE_SETTINGS[self.section]

		for param in params:
			default, _ = params[param]
			widget = self.widgets[param]

			if isinstance(widget, QLineEdit):
				widget.setText(default)
			elif isinstance(widget, RNASuiteBrowseLineEdit):
				widget.set_text(default)

class RNASuiteRGeneralSettingPage(RNASuiteGlobalSettingPage):
	section = 'R'

	def __init__(self, parent=None):
		super().__init__(parent)
		self.settings = QSettings()
		self.layout = QVBoxLayout()
		self.setLayout(self.layout)

	@property
	def params(self):
		return [
			AttrDict(
				option = 'binary',
				type = QFile,
				label = "Change R Binary"
			),
			AttrDict(
				option = 'cran_mirror',
				type = str,
				label = "CRAN Mirror"
			),
			AttrDict(
				option = 'mirror_list',
				type = list,
				label = None
			),
			AttrDict(
				option = 'mirror_custom',
				type = bool,
				text = "Allow custom cran mirror",
				label = None
			),
			AttrDict(
				option = 'bioc_mirror',
				type = str,
				label = "Bioconductor Mirror"
			)
		]

	def register_events(self):
		self.widgets.cran_mirror.setReadOnly(True)
		self.widgets.mirror_list.addItem(QListWidgetItem("Loading cran mirrors..."))
		self.widgets.mirror_custom.stateChanged.connect(self.on_custom_canr_mirror)
		self.load_cran_mirrors()

	def on_custom_canr_mirror(self, flag):
		self.widgets.cran_mirror.setReadOnly(not flag)

	def update_cran_mirrors(self, rows):
		self.mirror_urls = {}
		self.widgets.mirror_list.clear()

		for i, row in enumerate(rows[1:-1]):
			cols = row.split('\t')
			item = QListWidgetItem("{} - {}".format(cols[1], cols[5]))
			self.mirror_urls[i] = cols[4]
			self.widgets.mirror_list.addItem(item)

	@Slot()
	def on_mirror_changed(self, index):
		url = self.mirror_urls.get(index, '')
		self.widgets.cran_mirror.setText(url)

	def load_cran_mirrors(self):
		self.widgets.mirror_list.currentRowChanged.connect(self.on_mirror_changed)
		thread = RNASuiteCranMirrorThread(self)
		thread.result.connect(self.update_cran_mirrors)
		thread.start()

class RNASuiteContrastVersusWidget(QWidget):
	contrast_changed = Signal(int)

	def __init__(self, parent=None):
		super().__init__(parent)
		self.setMinimumSize(QSize(100, 150))

		self.create_widgets()
		self.set_layouts()

	def create_widgets(self):
		self.contrast_tree = QTreeWidget(self)
		self.contrast_tree.setRootIsDecorated(False)
		self.contrast_tree.setHeaderLabels(['Treatment', 'Control'])
		self.contrast_tree.header().setSectionResizeMode(QHeaderView.Stretch)
		self.contrast_tree.setContextMenuPolicy(Qt.CustomContextMenu)
		self.contrast_tree.customContextMenuRequested.connect(self.show_context_menu)
		self.treatment_select = QComboBox(self)
		self.control_select = QComboBox(self)
		self.versus_label = QLabel('vs', self)
		self.add_button = QPushButton(self)
		self.add_button.setText('Add')
		self.add_button.clicked.connect(self.on_add_clicked)

	def set_layouts(self):
		top_layout = QHBoxLayout()
		top_layout.setContentsMargins(0, 0, 0, 0)
		top_layout.addWidget(self.treatment_select, 1)
		top_layout.addWidget(self.versus_label)
		top_layout.addWidget(self.control_select, 1)
		top_layout.addWidget(self.add_button)

		main_layout = QVBoxLayout()
		main_layout.setContentsMargins(0, 0, 0, 0)
		main_layout.addLayout(top_layout)
		main_layout.addWidget(self.contrast_tree, 1)

		self.setLayout(main_layout)

	@Slot()
	def show_context_menu(self, pos):
		item = self.contrast_tree.currentItem()

		del_act = QAction("Delete")
		del_act.setDisabled(not item)
		del_act.triggered.connect(self.on_delete_contrast)
		
		clr_act = QAction("Clear")
		#clr_act.setDisabled(not item)
		clr_act.triggered.connect(self.on_clear_contrasts)
		
		menu = QMenu(self.contrast_tree)
		menu.addAction(del_act)
		menu.addAction(clr_act)
		menu.exec(self.contrast_tree.mapToGlobal(pos))

	@Slot()
	def on_delete_contrast(self):
		item = self.contrast_tree.currentItem()

		if not item:
			return

		index = self.contrast_tree.indexOfTopLevelItem(item)
		self.contrast_tree.takeTopLevelItem(index)
		self.contrast_changed.emit(self.contrast_tree.topLevelItemCount())

	@Slot()
	def on_clear_contrasts(self):
		self.contrast_tree.clear()
		self.contrast_changed.emit(0)

	@Slot()
	def on_add_clicked(self):
		treatment = self.treatment_select.currentText()
		control = self.control_select.currentText()

		if treatment and control:
			self.add_contrast(treatment, control)

	def add_contrast(self, treatment, control):
		item = QTreeWidgetItem(self.contrast_tree)
		item.setText(0, treatment)
		item.setText(1, control)
		self.contrast_tree.addTopLevelItem(item)
		self.contrast_changed.emit(self.contrast_tree.topLevelItemCount())

	def set_contrasts(self, contrasts):
		for treatment, control in contrasts:
			self.add_contrast(treatment, control)

	def get_contrasts(self):
		contrasts = []
		for i in range(self.contrast_tree.topLevelItemCount()):
			item = self.contrast_tree.topLevelItem(i)
			contrasts.append([item.text(0), item.text(1)])

		return contrasts

	def set_selection(self, items):
		self.treatment_select.addItems(items)
		self.control_select.addItems(items)

class RNASuiteColorGroups(QWidget):
	def __init__(self, parent=None, count=0):
		super().__init__(parent)
		self.color_counts = count
		self.color_widgets = []
		self.layout = QHBoxLayout()
		self.layout.setContentsMargins(0, 0, 0, 0)
		self.layout.addStretch()
		self.setLayout(self.layout)

		self.create_color_buttons()

	def add_color_button(self):
		widget = RNASuiteColorButton(self)
		self.color_widgets.append(widget)
		self.layout.addWidget(widget, 0, Qt.AlignLeft)
		return widget

	def remove_color_button(self):
		widget = self.color_widgets.pop()
		self.layout.removeWidget(widget)
		widget.deleteLater()

	def create_color_buttons(self):
		for i in range(self.color_counts):
			self.add_color_button()

	def change_color_buttons(self, count):
		if count > self.color_counts:
			d = count - self.color_counts

			for i in range(d):
				self.add_color_button()

		elif count < self.color_counts:
			d = self.color_counts - count

			for i in range(d):
				self.remove_color_button()

		self.color_counts = count

	def clear_colors(self):
		for i in range(self.color_counts):
			self.remove_color_button()

		self.color_counts = 0

	def set_colors(self, colors):
		self.clear_colors()

		for color in colors:
			widget = self.add_color_button()
			widget.set_color(color)

		self.color_counts = len(colors)

	def get_colors(self):
		return [w.get_color() for w in self.color_widgets]

