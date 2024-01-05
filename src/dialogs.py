from packaging.version import Version

from rchitect.utils import Rhome

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

from utils import *
from params import *
from widgets import *
from threads import *

__all__ = ['RNASuitePackageManagerDialog', 'RNASuiteDeseqParameterDialog',
	'RNASuiteShowDEGParameterDialog', 'RNASuiteGlobalSettingDialog',
	'RNASuiteEdgerParameterDialog', 'RNASuiteDEGDistPlotParameterDialog'
]

class RNASuiteGlobalSettingDialog(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.setWindowTitle("Global Settings")
		self.list_widget = QListWidget(self)
		self.list_widget.setIconSize(QSize(32, 32))
		self.stack_widget = QStackedWidget(self)

		self.list_widget.currentRowChanged.connect(self.stack_widget.setCurrentIndex)

		self.button_box = QDialogButtonBox(
			QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel
		)

		self.button_box.accepted.connect(self.save_settings)
		self.button_box.rejected.connect(self.reject)
		self.button_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self.reset_settings)

		main_layout = QVBoxLayout()
		main_layout.setSpacing(20)
		widget_layout = QHBoxLayout()
		widget_layout.setSpacing(20)
		widget_layout.addWidget(self.list_widget)
		widget_layout.addWidget(self.stack_widget, 1)
		main_layout.addLayout(widget_layout)
		main_layout.addWidget(self.button_box)
		self.setLayout(main_layout)

		self.create_pages()

	def sizeHint(self):
		return QSize(700, 400)

	def create_pages(self):
		pages = [
			('R General', 'icons/rlogo.svg', RNASuiteRGeneralSettingPage)
		]

		list_width = 0
		for text, icon, pager in pages:
			item = QListWidgetItem(QIcon(icon), text)
			self.list_widget.addItem(item)
			page = pager(self)
			self.stack_widget.addWidget(page)
			item_width = self.list_widget.visualItemRect(item).width()

			if item_width > list_width:
				list_width = item_width
		
		list_width += 10
		self.list_widget.setFixedWidth(list_width)

	def save_settings(self):
		for i in range(self.stack_widget.count()):
			widget = self.stack_widget.widget(i)
			widget.write_settings()

		self.accept()

	def reset_settings(self):
		for i in range(self.stack_widget.count()):
			widget = self.stack_widget.widget(i)
			widget.restore_settings()

class RNASuitePackageManagerDialog(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)

		self.setWindowTitle("R Package Manager")

		self.progress = RNASuiteWaitingSpinner(self)
		self.progress.hide()
		spacer = RNASuiteSpacerWidget(self)
		self.update_btn = QPushButton("Update", self)
		self.update_btn.setIcon(QIcon('icons/update.svg'))
		self.package_manager = RNASuitePackageTreeView(self)
		self.update_btn.clicked.connect(self.on_update_status)
		self.status_info = RNASuitePackageInstallMessage(self)
		self.package_manager.error.connect(self.on_error_occurred)
		self.package_manager.message.connect(self.on_update_message)
		self.package_manager.started.connect(self.on_install_started)
		self.package_manager.stopped.connect(self.on_install_stopped)

		self.layout = QVBoxLayout()
		self.setLayout(self.layout)
		top_layout = QHBoxLayout()
		top_layout.addWidget(self.progress)
		top_layout.addWidget(spacer)
		top_layout.addWidget(self.update_btn)
		self.layout.addLayout(top_layout)
		self.layout.addWidget(self.package_manager)
		self.layout.addWidget(self.status_info)

	def closeEvent(self, event):
		if self.package_manager.task_running():
			info = (
				"Are you sure you want to close the package manager?\n"
				"Closing the manager will stop the current package installation"
			)
			btn = QMessageBox.question(self, "Warnning", info)

			if btn == QMessageBox.Yes:
				self.package_manager.stop_task()
				event.accept()
			else:
				event.ignore()
				return

		event.accept()

	@Slot()
	def on_error_occurred(self, error):
		QMessageBox.critical(self, "Error", error)

	@Slot()
	def on_update_message(self, text):
		self.status_info.insertPlainText(text)
		scroll_bar = self.status_info.verticalScrollBar()
		scroll_bar.setValue(scroll_bar.maximum())

	@Slot()
	def on_update_status(self):
		if self.package_manager.task_running():
			return QMessageBox.warning(self, "Warnning", "A package installation is running")

		self.package_manager.update_version()

	@Slot()
	def on_install_started(self):
		self.progress.show()
		self.progress.start()

	@Slot()
	def on_install_stopped(self):
		self.progress.hide()
		self.progress.stop()

class RNASuiteParameterDialog(QDialog):
	parameters = {}
	title = None

	def __init__(self, parent=None, defines={}, dataset=None):
		super().__init__(parent)
		self.setWindowTitle(self.title)

		self.widgets = AttrDict()
		self.defines = defines
		self.dataset = dataset

		main_layout = QVBoxLayout()
		self.setLayout(main_layout)

		self.widget_layout = QFormLayout()

		self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		#self.buttons.button(QDialogButtonBox.Ok).setText("Start")
		self.buttons.accepted.connect(self.accept)
		self.buttons.rejected.connect(self.reject)

		main_layout.addLayout(self.widget_layout)
		main_layout.addWidget(self.buttons)

		self.register_widgets()
		self.register_events()

	def sizeHint(self):
		return QSize(400, 0)

	def register_widgets(self):
		for i, p in enumerate(self.parameters):
			val = self.defines.get(p.key, None) or p.default

			if p.type == int:
				self.widgets[p.key] = QSpinBox(self)
				self.widgets[p.key].setRange(*p.range)
				self.widgets[p.key].setSingleStep(p.step)
				self.widgets[p.key].setValue(val)

			elif p.type == float:
				self.widgets[p.key] = QDoubleSpinBox(self)
				self.widgets[p.key].setRange(*p.range)
				self.widgets[p.key].setSingleStep(p.step)
				self.widgets[p.key].setDecimals(5)
				self.widgets[p.key].setValue(val)

			elif p.type == str:
				self.widgets[p.key] = QLineEdit(self)

				if val:
					self.widgets[p.key].setText(val)

			elif p.type == list:
				self.widgets[p.key] = QComboBox(self)
				self.widgets[p.key].addItems(p.options)

				if isinstance(val, int):
					self.widgets[p.key].setCurrentIndex(val)

				elif val:
					self.widgets[p.key].setCurrentText(val)

			elif p.type == bool:
				self.widgets[p.key] = QCheckBox(self)

				if val:
					self.widgets[p.key].setCheckState(Qt.Checked)
				else:
					self.widgets[p.key].setCheckState(Qt.Unchecked)

			elif p.type == set:
				self.widgets[p.key] = RNASuiteMultipleSelect(self)
				self.widgets[p.key].add_items(p.options)

				if val:
					self.widgets[p.key].set_text(val)

			elif p.type == 'text':
				self.widgets[p.key] = QPlainTextEdit(self)

				if val:
					self.widgets[p.key].appendPlainText(val)

			elif p.type == 'color':
				self.widgets[p.key] = RNASuiteColorButton(self)

				if val:
					self.widgets[p.key].set_color(val)

			#label = QLabel(p.display, self)
			self.widget_layout.addRow(p.display, self.widgets[p.key])

			if 'help' in p:
				info = QLabel("<font color='gray'>{}</font>".format(p.help), self)
				#info.setWordWrap(True)
				self.widget_layout.addRow('', info)

	def register_events(self):
		pass

	def get_param_values(self):
		params = {}

		for k, w in self.widgets.items():
			if isinstance(w, QAbstractSpinBox):
				params[k] = w.value()

			elif isinstance(w, QLineEdit):
				params[k] = w.text().strip()

			elif isinstance(w, QComboBox):
				p = self.parameters[k]

				if p.index:
					params[k] = w.currentIndex()
				else:
					params[k] = w.currentText()

			elif isinstance(w, RNASuiteMultipleSelect):
				params[k] = w.get_text()

			elif isinstance(w, QCheckBox):
				params[k] = w.checkState() == Qt.Checked

		return params

	@classmethod
	def get_params(cls, parent=None, defines={}, dataset=None):
		dlg = cls(parent, defines, dataset)

		if dlg.exec() == QDialog.Accepted:
			return dlg.get_param_values()

class RNASuiteDeseqParameterDialog(RNASuiteParameterDialog):
	parameters = RNASuiteDEGParameters
	title = "Identify DEGs by DESeq2"

	@Slot()
	def custom_design_toggle(self, state):
		if state:
			self.widgets.design.setReadOnly(False)
		else:
			self.widgets.design.setReadOnly(True)

	@Slot()
	def comparison_between_changed(self, group):
		groups = list(map(str, self.dataset[group]))
		self.widgets.control.clear()
		self.widgets.control.addItems(groups)
		self.widgets.treatment.clear()
		self.widgets.treatment.addItems(groups)

		factors = list(self.dataset.keys())
		factors.remove(group)
		self.widgets.eliminate.add_items(factors)

		self.widgets.design.setText("~ {}".format(group))

	@Slot()
	def eliminate_factor_changed(self, factors):
		group = self.widgets.compare.currentText()

		if factors:
			self.widgets.design.setText("~ {} + {}".format(factors, group))
		else:
			self.widgets.design.setText("~ {}".format(group))

	def register_events(self):
		self.widgets.design.setReadOnly(True)
		self.widgets.custom.stateChanged.connect(self.custom_design_toggle)
		self.widgets.compare.currentTextChanged.connect(self.comparison_between_changed)
		self.widgets.compare.addItems(list(self.dataset.keys()))
		self.widgets.eliminate.currentTextChanged.connect(self.eliminate_factor_changed)

class RNASuiteEdgerParameterDialog(RNASuiteParameterDialog):
	parameters = RNASuiteEdgerParameters
	title = "Identify DEGs by edgeR"

	@Slot()
	def custom_design_toggle(self, state):
		if state:
			self.widgets.design.setReadOnly(False)
		else:
			self.widgets.design.setReadOnly(True)

	@Slot()
	def comparison_between_changed(self, group):
		groups = list(map(str, self.dataset[group]))
		self.widgets.control.clear()
		self.widgets.control.addItems(groups)
		self.widgets.treatment.clear()
		self.widgets.treatment.addItems(groups)

		factors = list(self.dataset.keys())
		factors.remove(group)
		self.widgets.eliminate.add_items(factors)

		self.widgets.design.setText("~ 0 + {}".format(group))

	@Slot()
	def eliminate_factor_changed(self, factors):
		group = self.widgets.compare.currentText()

		if factors:
			self.widgets.design.setText("~ 0 + {} + {}".format(group, factors))
		else:
			self.widgets.design.setText("~ 0 + {}".format(group))

	def with_replicate_changed(self, index):
		if index == 0:
			self.widgets.bcv.setDisabled(True)
			self.widgets.method.model().item(0).setEnabled(True)
			self.widgets.method.setCurrentIndex(0)
		else:
			self.widgets.bcv.setDisabled(False)
			self.widgets.method.model().item(0).setEnabled(False)
			self.widgets.method.setCurrentIndex(2)

	def register_events(self):
		self.widgets.design.setReadOnly(True)
		self.widgets.bcv.setDisabled(True)
		self.widgets.custom.stateChanged.connect(self.custom_design_toggle)
		self.widgets.compare.currentTextChanged.connect(self.comparison_between_changed)
		self.widgets.compare.addItems(list(self.dataset.keys()))
		self.widgets.eliminate.currentTextChanged.connect(self.eliminate_factor_changed)
		self.widgets.replicate.currentIndexChanged.connect(self.with_replicate_changed)

class RNASuiteShowDEGParameterDialog(RNASuiteParameterDialog):
	parameters = RNASuiteShowDEGParameters
	title = "Show DEGs"

	def register_events(self):
		group = self.defines['compare']
		groups = list(map(str, self.dataset[group]))
		self.widgets.control.clear()
		self.widgets.control.addItems(groups)
		self.widgets.treatment.clear()
		self.widgets.treatment.addItems(groups)

class RNASuiteDEGDistPlotParameterDialog(RNASuiteParameterDialog):
	parameters = RNASuiteDEGDistPlotParameters
	title = 'DEG Distribution Plot'
