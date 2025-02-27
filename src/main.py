import os
import sys
import multiprocessing

import rchitect

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from config import *
from window import *

class RNASuiteApplication(QApplication):
	osx_open_with = Signal(str)

	def __init__(self, argv):
		super().__init__(argv)
		self.set_style()

	def event(self, event):
		if sys.platform == 'darwin':
			if isinstance(event, QFileOpenEvent):
				self.osx_open_with.emit(event.file())

		return super().event(event)

	def set_style(self):
		with open('geneious.qss', encoding='utf-8') as style:
			sheet = style.read()

		self.setStyleSheet(sheet)

if __name__ == '__main__':
	multiprocessing.freeze_support()

	if os.name == 'nt':
		import ctypes
		myappid = "Dulab.RNASuite.RNASuite.{}".format(RNASUITE_VERSION)
		ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

	QCoreApplication.setOrganizationName("Dulab")
	QCoreApplication.setOrganizationDomain("big.cdu.edu.cn")
	QCoreApplication.setApplicationName("RNASuite")
	QCoreApplication.setApplicationVersion(RNASUITE_VERSION)
	QSettings.setDefaultFormat(QSettings.IniFormat)

	app = RNASuiteApplication(sys.argv)
	#r = RNASuiteREnvironment()
	win = RNASuiteMainWindow()

	#app.osx_open_with.connect(win.create_db_connect)

	args = app.arguments()

	if len(args) > 1:
		if os.path.isfile(args[v]):
			#win.open_project_file(args[1])
			pass

	sys.exit(app.exec())
	#rchitect.loop()
