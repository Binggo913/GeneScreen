from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon, QColor, QPainter, QPen, QPixmap
from PySide6.QtWidgets import QApplication, QWidget


def _device_pixel_ratio(widget: QWidget | None) -> float:
    if widget and widget.windowHandle():
        screen = widget.windowHandle().screen()
    else:
        screen = QApplication.primaryScreen()
    return screen.devicePixelRatio() if screen else 1.0


def create_icon_pixmap(size: int, widget: QWidget | None = None) -> QPixmap:
    dpr = _device_pixel_ratio(widget)
    pixmap = QPixmap(int(size * dpr), int(size * dpr))
    pixmap.setDevicePixelRatio(dpr)
    pixmap.fill(Qt.transparent)
    return pixmap


def draw_sidebar_icon(name: str, color: QColor, size: int = 20, widget: QWidget | None = None) -> QPixmap:
    pixmap = create_icon_pixmap(size, widget)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing)
    pen = QPen(color)
    pen.setWidthF(1.8)
    pen.setCapStyle(Qt.RoundCap)
    pen.setJoinStyle(Qt.RoundJoin)
    painter.setPen(pen)

    if name == "grid":
        margin = 3
        gap = 3
        cell = (size - margin * 2 - gap) / 2
        for row in range(2):
            for col in range(2):
                x = margin + col * (cell + gap)
                y = margin + row * (cell + gap)
                painter.drawRect(x, y, cell, cell)
    elif name == "id":
        margin = max(2.0, size * 0.15)
        width = size - margin * 2
        height = size - margin * 2
        radius = max(2.2, size * 0.18)
        painter.drawRoundedRect(margin, margin + 0.6, width, height - 1.2, radius, radius)
        dot = max(3.0, size * 0.22)
        painter.drawEllipse(margin + dot * 0.6, margin + dot * 0.8, dot, dot)
        line_x = margin + dot * 2.2
        line_y1 = margin + dot * 1.0
        line_y2 = margin + dot * 2.0
        painter.drawLine(line_x, line_y1, size - margin - dot * 0.6, line_y1)
        painter.drawLine(line_x, line_y2, size - margin - dot * 0.9, line_y2)
    elif name == "pin":
        center_x = size / 2
        center_y = size / 2 - 2
        radius = 4.5
        painter.drawEllipse(center_x - radius, center_y - radius, radius * 2, radius * 2)
        painter.drawLine(center_x, center_y + radius - 0.5, center_x, size - 3)
    elif name == "code":
        left = 5
        right = size - 5
        top = 5
        bottom = size - 5
        mid = size / 2
        painter.drawLine(left + 3, top, left, mid)
        painter.drawLine(left, mid, left + 3, bottom)
        painter.drawLine(right - 3, top, right, mid)
        painter.drawLine(right, mid, right - 3, bottom)
    elif name == "clock":
        radius = 7
        center = size / 2
        painter.drawEllipse(center - radius, center - radius, radius * 2, radius * 2)
        painter.drawLine(center, center, center, center - 4)
        painter.drawLine(center, center, center + 3.5, center + 2)
    elif name == "database":
        margin = 3
        width = size - margin * 2
        top_h = 4
        painter.drawEllipse(margin, margin, width, top_h * 2)
        painter.drawLine(margin, margin + top_h, margin, size - margin - top_h)
        painter.drawLine(margin + width, margin + top_h, margin + width, size - margin - top_h)
        painter.drawArc(margin, size - margin - top_h * 2, width, top_h * 2, 0, -180 * 16)
    else:
        y_positions = (5, 10, 15)
        for i, y in enumerate(y_positions):
            painter.drawLine(4, y, size - 4, y)
            knob_x = 7 if i % 2 == 0 else size - 7
            painter.drawEllipse(knob_x - 2.2, y - 2.2, 4.4, 4.4)

    painter.end()
    return pixmap


def build_sidebar_icon(
    name: str,
    normal_color: QColor,
    active_color: QColor,
    size: int = 20,
    widget: QWidget | None = None
) -> QIcon:
    normal = draw_sidebar_icon(name, normal_color, size=size, widget=widget)
    active = draw_sidebar_icon(name, active_color, size=size, widget=widget)
    icon = QIcon()
    icon.addPixmap(normal, QIcon.Normal, QIcon.Off)
    icon.addPixmap(active, QIcon.Normal, QIcon.On)
    return icon
