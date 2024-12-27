import os
import yaml

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from widgets import *

__all__ = [
	'RNASuitePluginManager'
	'RNASuitePluginRegister',
	'RNASuitePluginStarter',
	'RNASuitePluginWorker'
]

class RNASuitePluginManager(QObject):
	pass

class RNASuitePluginRegister(QObject):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.plugin_root = os.path.join(os.path.dirname(__file__), 'plugins')

	def run(self):
		it = QDirIterator(self.plugin_root, ['config.yaml'],
			QDir.Files, QDirIterator.Subdirectories)
		while it.hasNext():
			print(it.next())
		

class RNASuitePluginStarter(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)

class RNASuitePluginWorker(QObject):
	def __init__(self, parent=None):
		super().__init__(parent)

if __name__ == '__main__':
	RNASuitePluginRegister()
