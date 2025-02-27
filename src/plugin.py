import os
import yaml

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from models import *
from config import *
from widgets import *

__all__ = [
	'RNASuitePluginManager',
	'RNASuitePluginRegister',
	'RNASuitePluginStarter',
	'RNASuitePluginWorker'
]

class RNASuitePluginManager(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Plugins Manager")

		self.create_widgets()
		self.read_plugins()

	def sizeHint(self):
		return QSize(500, 300)

	def create_widgets(self):
		self.tree = QTreeView(self)
		self.tree.header().setStretchLastSection(False)
		self.model = RNASuitePluginModel(self)
		self.tree.setModel(self.model)

		button_box = QDialogButtonBox(QDialogButtonBox.Ok)
		button_box.accepted.connect(self.accept)

		main_layout = QVBoxLayout()
		main_layout.addWidget(self.tree)
		main_layout.addWidget(button_box)
		self.setLayout(main_layout)

	def update_column_width(self):
		self.tree.header().setSectionResizeMode(0, QHeaderView.ResizeToContents)
		self.tree.header().setSectionResizeMode(1, QHeaderView.Stretch)
		self.tree.header().setSectionResizeMode(2, QHeaderView.ResizeToContents)
		self.tree.header().setSectionResizeMode(3, QHeaderView.Interactive)


	def read_plugins(self):
		yaml_files = QDirIterator(RNASUITE_PLUGIN, ['config.yaml'],
			QDir.Files, QDirIterator.Subdirectories)

		while yaml_files.hasNext():
			yaml_file = yaml_files.next()

			with open(yaml_file) as stream:
				try:
					config = yaml.safe_load(stream)
				except yaml.YAMLError as err:
					print(err)
					continue

				self.model.add_plugin(config['name'], config['description'], config['version'], '')

		self.update_column_width()



class RNASuitePluginRegister(QObject):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.run()

	def run(self):
		it = QDirIterator(RNASUITE_PLUGIN, ['config.yaml'],
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
