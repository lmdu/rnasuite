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

__all__ = ['create_parameter_widget', 'get_widgets_parameters',
	'RNASuiteInputListWidget', 'RNASuiteOutputTreeWidget',
	'RNASuitePackageInstallButton', 'RNASuiteWaitingSpinner',
	'RNASuitePackageTreeView', 'RNASuitePackageInstallMessage',
	'RNASuiteSpacerWidget', 'RNASuiteMultipleSelect',
	'RNASuiteRGeneralSettingPage', 'RNASuiteColorButton',
	'RNASuiteContrastVersusWidget', 'RNASuiteColorGroups'
]

def create_parameter_widget(param, value=None):
	match param.type:
		case 'int':
			widget = QSpinBox()
			widget.setRange(*param.range)
			widget.setSingleStep(param.step)

			if value is not None:
				widget.setValue(value)

		case 'float':
			widget = QDoubleSpinBox()
			widget.setRange(*param.range)
			widget.setSingleStep(param.step)
			widget.setDecimals(5)

			if value is not None:
				widget.setValue(value)

		case 'str':
			widget = QLineEdit()

			if value is not None:
				widget.setText(value)

		case 'list':
			widget = QComboBox()
			widget.addItems(param.options)

			if value is not None:
				if isinstance(value, int):
					widget.setCurrentIndex(value)

				else:
					widget.setCurrentText(value)

		case 'bool':
			widget = QCheckBox()

			if value is not None:
				if value:
					widget.setCheckState(Qt.Checked)

				else:
					widget.setCheckState(Qt.Uncheched)

		case 'select':
			widget = RNASuiteMultipleSelect()
			widget.add_items(param.options)

			if value is not None:
				widget.set_text(value)

		case 'text':
			widget = QPlainTextEdit()

			if value is not None:
				widget.appendPlainText(value)

		case 'color':
			widget = RNASuiteColorButton()

			if value is not None:
				widget.set_color(value)

		case 'colors':
			widget = RNASuiteColorGroups()

			if value is not None:
				widget.set_colors(value)

		case 'contrast':
			widget = RNASuiteContrastVersusWidget()

			if value is not None:
				widget.set_contrasts(value)

	return widget

def get_widgets_parameters(widgets, params):
	values = {}

	for k, w in widgets.items():
		p = params[k]

		if not p.expose:
			continue

		match w:
			case QAbstractSpinBox():
				values[k] = w.value()

			case QLineEdit():
				values[k] = w.text().strip()

			case QComboBox():
				if p.index:
					values[k] = w.currentIndex()
				else:
					values[k] = w.currentText()

			case QCheckBox():
				values[k] = w.checkState() == Qt.Checked

			case QPlainTextEdit():
				values[k] = w.toPlainText()

			case RNASuiteMultipleSelect():
				values[k] = w.get_text()

			case RNASuiteColorButton():
				values[k] = w.get_color()

			case RNASuiteColorGroups():
				values[k] = w.get_colors()

			case RNASuiteContrastVersusWidget():
				values[k] = w.get_contrasts()

	return values

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
		self.itemDoubleClicked.connect(self._on_open_table)

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
	show_table = Signal(pandas.DataFrame)
	show_panel = Signal(str)
	show_plot = Signal(int)

	def __init__(self, parent=None):
		super().__init__(parent)
		self.datasets = {}
		self.setRootIsDecorated(False)
		self.doubleClicked.connect(self._on_row_clicked)

		self.create_model()

	def sizeHint(self):
		return QSize(200, 500)

	def create_model(self):
		self._data = pandas.DataFrame(columns=['type', 'name', 'update', 'plot'])
		self._model = RNASuiteOutputTreeModel(self)
		self._model.load_data(self._data)
		self.setModel(self._model)
		self.setColumnHidden(3, True)
		self.header().setStretchLastSection(False)
		self.header().setSectionResizeMode(0, QHeaderView.ResizeToContents)
		self.header().setSectionResizeMode(1, QHeaderView.Stretch)
		self.header().setSectionResizeMode(2, QHeaderView.ResizeToContents)

	def add_row(self, **row):
		key = "{}-{}-{}".format(row['type'], row['name'], row['plot'])

		if key in self.datasets:
			index = self._data[(self._data['type'] == row['type']) & (self._data['name'] == row['name']) & (self._data['plot'] == row['plot'])].index
			self._data.iloc[index, 2] = 1

		else:
			self._data = pandas.concat(
				[pandas.DataFrame([row], columns = self._data.columns), self._data],
				ignore_index = True
			)
			self._model.load_data(self._data)
		
		if row['plot']:
			chart = row.pop('chart')
			self.datasets[key] = chart

		else:
			data = row.pop('data')
			data_frame = pandas.DataFrame.from_dict(data, orient='tight')
			self.datasets[key] = data_frame

	@Slot()
	def receive(self, rtype, result):
		match rtype:
			case 'degs':
				v = result['degs_versus']
				contrast = "{} vs {}".format(v[-2], v[-1])

				if 'normal_count' in result:
					self.add_row(
						type = 'Genes',
						name = 'Noramlized Counts',
						plot = 0,
						update = 1,
						data = result['normal_count']
					)

				if 'degs_list' in result:
					self.add_row(
						type = 'DEGs',
						name = contrast,
						plot = 0,
						update = 1,
						data = result['degs_list']
					)

				if 'degs_plot' in result:
					self.add_row(
						type = result['plot_type'],
						name = contrast,
						plot = 1,
						update = 0,
						chart = result['degs_plot']
					)

					self.show_panel.emit('deseq_maplot')

	def _on_row_clicked(self, index):
		row = index.row()
		_type = self._data.iloc[row, 0]
		name = self._data.iloc[row, 1]
		plot = self._data.iloc[row, 3]
		key = "{}-{}-{}".format(_type, name, plot)

		if plot:
			pass

		else:
			data = self.datasets[key]
			self.show_table.emit(data)
			self._data.iloc[row, 2] = 0
			#self._model.refresh()

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

	def set_colors(self, colors):
		for color in colors:
			widget = self.add_color_button()
			widget.set_color(color)

		self.color_counts = len(colors)

	def get_colors(self):
		return [w.get_color() for w in self.color_widgets]

