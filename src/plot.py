import json
import time

from PySide6.QtSvg import *
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

__all__ = [
	'RNASuitePlotViewer',
]

class RNASuitePlotControlPanel(QWidget):
	parameters = None
	function = None
	plotname = None

	def __init__(self, parent=None):
		super().__init__(parent)
		self.rid = None
		self.pyid = None
		self.parent = parent
		self.widgets = AttrDict()
		self.stored = {}

		self.create_widgets()
		self.set_layouts()
		self.register_widgets()
		self.register_events()
		self.group_widgets()

	def set_plot(self, rid, pyid):
		self.rid = rid
		self.pyid = pyid

		params = self.stored.get(self.pyid, None)

		if params is not None:
			self.reset_widgets(params)

		else:
			self.restore_widgets()

		#self.title_label.setText("<b>{} {}</b>".format(self.plotname, self.rid))

	@Slot()
	def _on_update_clicked(self):
		if self.parent.has_running_worker():
			return

		params = self.get_param_values()
		self.stored[self.pyid] = params
		params['id'] = self.rid

		self.parent.pyconn.send({
			'action': 'call',
			'func': self.function,
			'params': params
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
		self.main_layout.setContentsMargins(3, 5, 0, 5)
		self.main_layout.addWidget(self.title_label)
		self.main_layout.addWidget(self.param_widgets, 1)
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

				if p.type == 'bool':
					self.widgets[k].setText(p.display.rstrip(':'))
					label = None
				else:
					label = QLabel(p.display)

				widgets.append((label, self.widgets[k]))

			self.param_widgets.add_accordions(group, widgets)

		self.param_widgets.add_stretcher()

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

class RNASuiteDegsVennPlotControlPanel(RNASuitePlotControlPanel):
	parameters = RNASuiteDegsVennPlotControlParameters
	function = 'rnasuite_degs_venn_plot_update'
	plotname = 'DEGs Venn Plot'

	@property
	def groups(self):
		return {
			'Show percent': ['show_percentage', 'digits'],
			'Fill colors': ['fill_color', 'fill_alpha'],
			'Stroke style': ['stroke_color', 'stroke_alpha', 'stroke_size', 'stroke_linetype'],
			'Set name style': ['set_name_color', 'set_name_size'],
			'Text style': ['text_color', 'text_size']
		}

class RNASuiteDegsUpsetPlotControlPanel(RNASuitePlotControlPanel):
	parameters = RNASuiteDegsUpsetPlotControlParameters
	function = 'rnasuite_degs_upset_plot_update'
	plotname = 'DEGs Upset Plot'

	@property
	def groups(self):
		return {
			'Main bar plot': ['main_ylabel', 'main_color', 'show_number', 'number_angles'],
			'Set bar plot': ['set_xlabel', 'set_color'],
			'Matrix plot': ['point_size', 'line_width', 'matrix_color', 'matrix_alpha'],
			'Order by': ['order_by', 'order_desc'],
			'Intersection': ['intersect_num', 'show_empty'],
			'Text size': ['text_scale']
		}

	def correct_param_values(self, values):
		if values['show_empty'] == 'off':
			values['show_empty'] = None

		return values

class RNASuiteDegsVolcanoPlotControlPanel(RNASuitePlotControlPanel):
	parameters = RNASuiteDegsVolcanoPlotControlParameters
	function = 'rnasuite_degs_volcano_plot_update'
	plotname = 'DEGs Volcano Plot'

	@property
	def groups(self):
		return {
			'Show labels': ['top'],
			'Point style': ['point_color', 'point_size'],
			'Show vertical line': ['show_vline', 'vline_type', 'vline_color', 'vline_width'],
			'Show horizontal line': ['show_hline', 'hline_type', 'hline_color', 'hline_width'],
			'Legend': ['legend_position'],
			'X limit': ['x_limit'],
			'Theme': ['theme_name', 'base_size']
		}

class RNASuitePlotController(QStackedWidget):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent
		self.panel_mapping = {}

		self.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)

	def sizeHint(self):
		return QSize(200, 10)

	def add_panel(self, panel):
		match panel:
			case 'deseq_maplot':
				panel_widget = RNASuiteDeseqMaPlotControlPanel(self.parent)

			case 'deg_distplot':
				panel_widget = RNASuiteDegsDistPlotControlPanel(self.parent)

			case 'deg_vennplot':
				panel_widget = RNASuiteDegsVennPlotControlPanel(self.parent)

			case 'deg_volcanoplot':
				panel_widget = RNASuiteDegsVolcanoPlotControlPanel(self.parent)

			case 'deg_upsetplot':
				panel_widget = RNASuiteDegsUpsetPlotControlPanel(self.parent)

		index = self.addWidget(panel_widget)
		self.panel_mapping[panel] = index

		return index

	def show_panel(self, panel, rid, pyid):
		if panel in self.panel_mapping:
			index = self.panel_mapping[panel]

		else:
			index = self.add_panel(panel)

		self.setCurrentIndex(index)
		self.widget(index).set_plot(rid, pyid)

class RNASuitePlotViewer(QWidget):
	error = Signal(str)

	def __init__(self, parent=None):
		super().__init__(parent)
		self.parent = parent

		self.plot_id = -1
		self.plot_num = 0
		self.plot_zoom = 1.0
		self.plot_width = 800
		self.plot_height = 600
		self.aspect_ratio = 0

		self.tool_bar = QToolBar(self)
		self.tool_bar.setIconSize(QSize(16, 16))
		
		self.zoom_label = QLabel("{:.1f}".format(self.plot_zoom), self)
		
		self.width_spin = QSpinBox(self)
		self.width_spin.setRange(1, 10000)
		self.width_spin.setValue(800)
		self.width_spin.setSingleStep(10)
		self.width_spin.valueChanged.connect(self.on_width_changed)
		self.width_spin.editingFinished.connect(self.redraw_plot)
		
		self.height_spin = QSpinBox(self)
		self.height_spin.setRange(1, 10000)
		self.height_spin.setValue(600)
		self.height_spin.setSingleStep(10)
		self.height_spin.valueChanged.connect(self.on_height_changed)
		self.height_spin.editingFinished.connect(self.redraw_plot)
		
		#self.ratio_check = QPushButton(self)
		#self.ratio_check.setCheckable(True)
		#self.ratio_check.setIcon(QIcon('icons/link.svg'))
		#self.ratio_check.toggled.connect(self.on_aspect_ratio_changed)

		#self.view = QSvgWidget(self)
		#self.plot_render = QSvgRenderer('icons/logo.svg')
		#self.plot_view = QGraphicsView()
		#self.view.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		#self.view.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		#self.view.setDragMode(QGraphicsView.ScrollHandDrag)
		#self.plot_scene = QGraphicsScene()
		#self.plot_view.setScene(self.plot_scene)

		#self.plot_item = QGraphicsSvgItem()
		#self.plot_item.setFlags(QGraphicsItem.ItemIsSelectable | QGraphicsItem.ItemIsMovable)
		#self.plot_item.setSharedRenderer(self.plot_render)
		#self.plot_scene.addItem(self.plot_item)
		self.plot_render = QSvgWidget(self)
		self.plot_view = QScrollArea(self)
		self.plot_view.setAlignment(Qt.AlignCenter)
		self.plot_view.setWidgetResizable(True)
		#self.plot_view.setBackgroundRole(QPalette.Light)
		self.plot_view.setWidget(self.plot_render)
		palette = self.plot_view.palette()
		palette.setColor(QPalette.Window, Qt.white)
		self.plot_view.setPalette(palette)

		self.plot_control = RNASuitePlotController(self)
		self.show_panel = self.plot_control.show_panel

		self.plot_split = QSplitter(self)
		self.plot_split.addWidget(self.plot_view)
		self.plot_split.addWidget(self.plot_control)
		self.plot_split.setStretchFactor(0, 1)

		layout = QVBoxLayout()
		layout.setContentsMargins(0, 0, 0, 0)
		layout.addWidget(self.tool_bar)
		layout.addWidget(self.plot_split)
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
			QLabel("Scale:", self),
			self.zoom_label,
			('icons/zoomin.svg', 'Zoom in', self.zoom_plot_in, False),
			('icons/zoomout.svg', 'Zoom out', self.zoom_plot_out, False),
			None,
			QLabel("Width:", self),
			self.width_spin,
			('icons/link.svg', 'Keep aspect ratio', self.on_aspect_ratio_changed, True),
			QLabel("Height:", self),
			self.height_spin,
			None,
			('icons/save.svg', 'Save plot to file', self.save_plot_file, False)
		]

		spacer = QWidget()
		spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
		self.tool_bar.addWidget(spacer)

		for item in actions:
			if item is None:
				self.tool_bar.addSeparator()

			elif isinstance(item, QWidget):
				self.tool_bar.addWidget(item)

			else:
				icon, text, func, toggle = item
				action = QAction(QIcon(icon), text, self)
				action.triggered.connect(func)
				action.setCheckable(toggle)
				self.tool_bar.addAction(action)

	#def resizeEvent(self, event):
	#	self.redraw_plot(index=-1)

	@Slot()
	def on_aspect_ratio_changed(self, checked):
		if checked:
			self.aspect_ratio = self.width_spin.value()/self.height_spin.value()
		else:
			self.aspect_ratio = 0

	@Slot()
	def on_width_changed(self, width):
		if self.aspect_ratio > 0:
			height = int(width/self.aspect_ratio)

			with QSignalBlocker(self.height_spin):
				self.height_spin.setValue(height)

			self.plot_height = height

		self.plot_width = width

	@Slot()
	def on_height_changed(self, height):
		if self.aspect_ratio > 0:
			width = int(height * self.aspect_ratio)

			with QSignalBlocker(self.width_spin):
				self.width_spin.setValue(width)

			self.plot_width = width

		self.plot_height = height

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

	def get_plot(self, **kwargs):
		url = QUrl("http://localhost:{}/plot".format(RNASUITE_PORT))
		query = QUrlQuery()
		query.addQueryItem('token', RNASUITE_TOKEN)
		query.addQueryItem('renderer', 'svgp')
		query.addQueryItem('id', str(self.plot_id))
		query.addQueryItem('zoom', str(self.plot_zoom))
		query.addQueryItem('width', str(self.plot_width))
		query.addQueryItem('height', str(self.plot_height))

		for k, v in kwargs.items():
			query.addQueryItem(k, str(v))

		url.setQuery(query)
		request = QNetworkRequest()
		request.setUrl(url)
		self.manager.get(request)

	@Slot()
	def on_error_occurred(self, error):
		self.error.emit(str(error))

	@Slot()
	def show_plot(self, static_id):
		self.plot_id = static_id
		self.get_plot()

	@Slot()
	def redraw_plot(self):
		#self.get_plot(index = -1)
		self.get_plot()

	@Slot()
	def on_text_received(self, text):
		data = json.loads(text)
		self.plot_num = data['hsize']
		self.redraw_plot()

	@Slot()
	def update_image(self, svg):
		if not isinstance(svg, QByteArray):
			svg = QByteArray(svg)

		self.plot_render.load(svg)
		self.plot_render.setFixedSize(self.plot_render.renderer().defaultSize())
		#self.plot_item.setSharedRenderer(self.plot_render)
		#self.plot_item.update()


		#self.plot_item.setSharedRenderer(self.plot_render)
		#item = QGraphicsSvgItem()
		#item.setFlags(QGraphicsItem.ItemIsSelectable | QGraphicsItem.ItemIsMovable)
		#item.setSharedRenderer(self.render)
		#self.scene.clear()
		#self.scene.addItem(item)

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

