import os
import math

from rchitect.utils import Rhome

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *
from PySide6.QtSvgWidgets import *

from utils import *
from config import *
from threads import *

__all__ = ['RNASuitePackageInstallButton', 'RNASuiteWaitingSpinner',
	'RNASuitePackageTreeView', 'RNASuitePackageInstallMessage',
	'RNASuiteSpacerWidget', 'RNASuiteMultipleSelect',
	'RNASuiteRGeneralSettingPage'
]

class RNASuiteSpacerWidget(QWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

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
