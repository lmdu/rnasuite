import json
import time

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *
from PySide6.QtNetwork import *
from PySide6.QtSvgWidgets import *
from PySide6.QtWebSockets import *

from utils import *
from config import *
from params import *
from widgets import *

__all__ = ['RNASuitePlotViewer', 'RNASuitePlotStackedWidget']

class RNASuitePlotViewer(QWidget):
	error = Signal(str)

	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent

		self.plot_id = 0
		self.plot_num = 0
		self.plot_zoom = 1.0

		self.tool = QToolBar(self)
		self.tool.setIconSize(QSize(16, 16))
		self.plot_label = QLabel("{} / {}".format(self.plot_id, self.plot_num), self)
		self.zoom_label = QLabel("{:.1f}".format(self.plot_zoom), self)

		self.view = QSvgWidget(self)
		layout = QVBoxLayout()
		layout.setContentsMargins(0, 0, 0, 0)
		#layout.addWidget(self.tool)
		layout.addWidget(self.view)
		self.setLayout(layout)

		self.socket = QWebSocket()
		self.socket.error.connect(self.on_error_occurred)
		self.socket.connected.connect(self.on_connected)
		self.socket.disconnected.connect(self.on_disconnected)
		self.socket.textMessageReceived.connect(self.on_text_received)
		#self.connect()

		self.manager = QNetworkAccessManager(self)
		self.manager.finished.connect(self.on_request_finished)

		self.create_actions()

	def create_actions(self):
		actions = [
			self.plot_label,
			('icons/prev.svg', 'Show previous plot', self.display_prev_plot),
			('icons/next.svg', 'Show next plot', self.display_next_plot),
			None,
			self.zoom_label,
			('icons/zoomin.svg', 'Zoom in', self.zoom_plot_in),
			('icons/zoomout.svg', 'Zoom out', self.zoom_plot_out),
			None,
			('icons/delete.svg', 'Remove current plot', self.remove_current_plot),
			('icons/clear.svg', 'Remove all plots', self.remove_all_plots),
			None,
			('icons/save.svg', 'Save plot to file', self.save_plot_file)
		]

		for item in actions:
			if item is None:
				self.tool.addSeparator()

			elif isinstance(item, QWidget):
				self.tool.addWidget(item)

			else:
				icon, text, func = item
				action = QAction(QIcon(icon), text, self)
				action.triggered.connect(func)
				self.tool.addAction(action)

	def resizeEvent(self, event):
		self.redraw_plot()

	@Slot()
	def on_connected(self):
		pass

	@Slot()
	def on_disconnected(self):
		pass

	@Slot()
	def connect_to_socket(self):
		url = QUrl("ws://localhost:{}".format(RNASUITE_PORT))
		self.socket.open(url)

	def disconnect_to_socket(self):
		self.socket.close()

	def zoom_plot_in(self):
		if not self.plot_num:
			return

		self.plot_zoom += 0.2

		if self.plot_zoom > 10:
			self.plot_zoom = 10

		self.redraw_plot()

		self.zoom_label.setText("{:.1f}".format(self.plot_zoom))

	def zoom_plot_out(self):
		if not self.plot_num:
			return

		self.plot_zoom -= 0.2

		if self.plot_zoom < 0.1:
			self.plot_zoom = 0.2

		self.redraw_plot()
		self.zoom_label.setText("{:.1f}".format(self.plot_zoom))

	def remove_current_plot(self):
		pass

	def remove_all_plots(self):
		pass

	def display_prev_plot(self):
		self.plot_id -= 1

		if self.plot_id < 0:
			self.plot_id = 0

		self.redraw_plot()

	def display_next_plot(self):
		self.plot_id += 1

		if self.plot_id == self.plot_num:
			self.plot_id = self.plot_num - 1

		self.redraw_plot()

	def save_plot_file(self):
		pass

	def remove_plot(self, static_id):
		#return static_id
		url = QUrl("http://localhost:{}/remove".format(RNASUITE_PORT))
		query = QUrlQuery()
		query.addQueryItem('token', RNASUITE_TOKEN)
		query.addQueryItem('id', str(static_id))
		url.setQuery(query)
		request = QNetworkRequest()
		request.setUrl(url)
		self.manager.get(request)

		#self.parent.pyconn.send({
		#	'action': 'remove',
		#	'params': {'page': plot_id+1}
		#})

	def show_plot(self, static_id):
		self.redraw_plot(static_id=static_id)

	def get_plot(self, **kwargs):
		url = QUrl("http://localhost:{}/plot".format(RNASUITE_PORT))
		query = QUrlQuery()
		query.addQueryItem('token', RNASUITE_TOKEN)
		query.addQueryItem('renderer', 'svgp')

		for k, v in kwargs.items():
			query.addQueryItem(k, str(v))

		url.setQuery(query)
		request = QNetworkRequest()
		request.setUrl(url)
		self.manager.get(request)

		#kwargs['renderer'] = 'svgp'
		#self.parent.pyconn.send({
		#	'action': 'plot',
		#	'params': kwargs
		#})

	@Slot()
	def on_error_occurred(self, error):
		self.error.emit(str(error))

	@Slot()
	def redraw_plot(self, index=None, static_id=None):
		if not self.plot_num:
			return

		size = self.size()
		width = int(size.width() * 0.95)
		height = int(size.height() * 0.95)

		if static_id is not None:
			self.get_plot(
				id = static_id,
				#page = self.plot_id + 1,
				zoom = self.plot_zoom,
				width = width,
				height = height
			)

		elif index is not None:
			self.get_plot(
				index = index,
				zoom = self.plot_zoom,
				width = width,
				height = height
			)

		self.plot_label.setText("{} / {}".format(self.plot_id+1, self.plot_num))

	@Slot()
	def on_text_received(self, text):
		data = json.loads(text)

		#if data['hsize'] > self.plot_num:
		#	self.plot_num = data['hsize']
		#	last_id = self.plot_num - 1

		#	print(last_id)
		self.plot_num = data['hsize']
		self.redraw_plot(index=-1)

	@Slot()
	def update_image(self, svg):
		if isinstance(svg, QByteArray):
			self.view.load(svg)
		else:
			self.view.load(QByteArray(svg))

	@Slot()
	def on_request_finished(self, reply):
		if reply.url().path() == '/remove':
			return

		if reply.error() == QNetworkReply.NoError:
			data = reply.readAll()
			self.update_image(data)

		else:
			self.error.emit(reply.errorString())

		reply.deleteLater()

class RNASuitePlotControlPanel(QWidget):
	parameters = None
	function = None
	plotname = None

	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent
		self.widgets = AttrDict()

		self.create_widgets()
		self.set_layouts()
		self.register_widgets()
		self.register_events()
		self.group_widgets()

	@Slot()
	def _on_update_clicked(self):
		self.parent.pyconn.send({
			'action': 'call',
			'func': self.function,
			'params': self.get_param_values()
		})

	def create_widgets(self):
		self.title_label = QLabel("<b>{}</b>".format(self.plotname), self)
		self.param_widgets = RNASuiteAccordionWidget(self)
		self.update_button = QPushButton(self)
		self.update_button.setText('Update Plot')
		self.update_button.setIcon(QIcon('icons/update.svg'))
		self.update_button.clicked.connect(self._on_update_clicked)
		self.reset_button = QPushButton(self)
		self.reset_button.setText('Restore Defaults')
		self.reset_button.clicked.connect(self.restore_widgets)

	def set_layouts(self):
		self.main_layout = QVBoxLayout(self)
		self.main_layout.addWidget(self.title_label)
		self.main_layout.addWidget(self.param_widgets)
		#self.main_layout.addLayout(self.widget_layout)
		self.main_layout.addWidget(self.update_button)
		self.main_layout.addWidget(self.reset_button)
		self.main_layout.addStretch()

	def register_widgets(self):
		for p in self.parameters:
			self.widgets[p.key] = create_parameter_widget(p)
			set_parameter_widget_value(self.widgets[p.key], p.default, p.index)

	def restore_widgets(self):
		for p in self.parameters:
			set_parameter_widget_value(self.widgets[p.key], p.default, p.index)

	def reset_widgets(self, params):
		for k, v in params.items():
			if k in self.widgets:
				p = self.parameters[k]
				set_parameter_widget_value(self.widgets[k], v, p.index)

	def group_widgets(self):
		for group in self.groups:
			keys = self.groups[group]

			widgets = []
			for k in keys:
				p = self.parameters[k]
				label = QLabel(p.display)
				widgets.append((label, self.widgets[k]))

			self.param_widgets.add_accordions(group, widgets)

	def register_events(self):
		pass

	def correct_param_values(self, values):
		return values

	def get_param_values(self):
		values = {}

		for k, w in self.widgets.items():
			p = self.parameters[k]

			if not p.expose:
				continue

			values[k] = get_parameter_widget_value(w, p.index)

		return self.correct_param_values(values)

class RNASuiteDeseqMaPlotControlPanel(RNASuitePlotControlPanel):
	parameters = RNASuiteDeseqMaPlotControlParameters
	function = 'rnasuite_deseq_ma_plot_update'
	plotname = 'DESeq2 MA Plot'

	@property
	def groups(self):
		return {
			'Title and labels': ['main', 'xlab'],
			'Y limit': ['ylim'],
			'Point colors': ['colNonSig', 'colSig'],
			'Line color': ['colLine']
		}

	def correct_param_values(self, values):
		if not any(values['ylim']):
			values['ylim'] = None

		degparam = self.parent.global_params['degs']
		values['treatment'] = degparam['treatment']
		values['control'] = degparam['control']

		return values

class RNASuiteDegsDistPlotControlPanel(RNASuitePlotControlPanel):
	parameters = RNASuiteDegsDistPlotControlParameters
	function = 'rnasuite_degs_dist_plot_update'
	plotname = 'DEGs Distribution Plot'

	@property
	def groups(self):
		return {
			'Plot type': ['plot_type'],
			'Title and labels': ['plot_title', 'x_label', 'x_rotate', 'y_label', 'show_label'],
			'Legend': ['legend_title'],
			'Fill colors': ['bar_colors'],
			'Theme': ['theme_name', 'base_size']
		}

class RNASuitePlotStackedWidget(QStackedWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent
		self.panel_mapping = {}

	def sizeHint(self):
		return QSize(200, 10)

	def add_panel(self, panel):
		match panel:
			case 'deseq_maplot':
				panel_widget = RNASuiteDeseqMaPlotControlPanel(self.parent)

			case 'deg_distplot':
				panel_widget = RNASuiteDegsDistPlotControlPanel(self.parent)

		index = self.addWidget(panel_widget)
		self.panel_mapping[panel] = index

		return index

	def show_panel(self, panel):
		if panel in self.panel_mapping:
			index = self.panel_mapping[panel]

		else:
			index = self.add_panel(panel)

		self.setCurrentIndex(index)
