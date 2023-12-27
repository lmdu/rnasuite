import os
import math

from rchitect.utils import Rhome

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *
from PySide6.QtSvgWidgets import *

from utils import *
from config import *

__all__ = ['RNASuitePackageInstallButton', 'RNASuiteWaitingSpinner',
	'RNASuitePackageTreeView', 'RNASuitePackageInstallMessage',
	'RNASuiteSpacerWidget', 'RNASuiteMultipleSelect'
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
		self.get_versions()
		self.create_model()

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

	def get_versions(self):
		self.packages = RNASuitePackages()
		rcode = (
			"rv<-list(name='R', ver=packageVersion('base'));"
			"write.table(rv, sep='\t', quote=F, row.names=F, col.names=F);"
			"ps<-installed.packages();"
			"ps<-ps[,c('Package','Version')];"
			"write.table(ps, sep='\t', quote=F, row.names=F, col.names=F)"
		)
		proc = QProcess()
		proc.errorOccurred.connect(self.on_error_occurred)
		proc.start(self.rscript, ['-e', rcode])
		proc.waitForFinished()
		res = proc.readAllStandardOutput()
		lines = res.data().decode()
		versions = {}

		for line in lines.split('\n'):
			if line.strip():
				p, v = line.split()
				versions[p.lower()] = v

		for c in self.packages.get_order():
			for p in self.packages[c]:
				p.version = versions.get(p.name.lower(), '')
				p.status = compare_version(p.required, p.version)

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
		self.update_data()

	def update_data(self):
		self._model.setHorizontalHeaderLabels([
			'Name', 'Description', 'Required', 'Version', 'Status'
		])
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

		rcode = None
		if package == 'R':
			QDesktopServices.openUrl("https://www.r-project.org/")
		elif repository == 'CRAN':
			rcode = "install.packages('{}')".format(package)
		elif repository == 'Bioconductor':
			rcode = "if (!require('BiocManager',quietly=TRUE)){{install.packages('BiocManager')}}; BiocManager::install('{}');".format(package)

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
		self.get_versions()
		self._model.clear()
		self.update_data()

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

