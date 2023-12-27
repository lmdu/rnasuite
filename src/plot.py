import json
import time

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *
from PySide6.QtNetwork import *
from PySide6.QtSvgWidgets import *
from PySide6.QtWebSockets import *

from config import *

__all__ = ['RNASuitePlotViewer']

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
		layout.addWidget(self.tool)
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
		self.url = QUrl("http://localhost:{}/plot".format(RNASUITE_PORT))

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

	def get_plot(self, **kwargs):
		#kwargs['token'] = HTTPGD_TOKEN
		kwargs['renderer'] = 'svgp'
		#query = '&'.join('{}={}'.format(k, v) for k, v in kwargs.items())
		#self.url.setQuery(query)
		#request = QNetworkRequest()
		#request.setUrl(self.url)
		#self.manager.get(request)
		self.parent.pyconn.send({
			'action': 'plot',
			'params': kwargs
		})

	@Slot()
	def on_error_occurred(self, error):
		self.error.emit(str(error))

	@Slot()
	def redraw_plot(self):
		if not self.plot_num:
			return

		size = self.size()
		width = int(size.width() * 0.95)
		height = int(size.height() * 0.95)

		self.get_plot(
			#id = self.plot_id,
			page = self.plot_id + 1,
			zoom = self.plot_zoom,
			width = width,
			height = height
		)

		self.plot_label.setText("{} / {}".format(self.plot_id+1, self.plot_num))

	@Slot()
	def on_text_received(self, text):
		data = json.loads(text)

		if data['hsize'] > self.plot_num:
			self.plot_num = data['hsize']
			self.plot_id = self.plot_num - 1
			self.redraw_plot()

	@Slot()
	def update_image(self, svg):
		if isinstance(svg, QByteArray):
			self.view.load(svg)
		else:
			self.view.load(QByteArray(svg))

	@Slot()
	def on_request_finished(self, reply):
		if reply.error() == QNetworkReply.NoError:
			data = reply.readAll()
			self.update_image(data)

		else:
			self.error.emit(reply.errorString())

		reply.deleteLater()

