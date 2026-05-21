"""Shared table checkbox helpers for selectable table columns."""

from typing import Optional

from PySide6.QtCore import QPoint, QPointF, Qt, QRect, QRectF, QSize
from PySide6.QtGui import QPainter, QPen
from PySide6.QtWidgets import (
    QCheckBox,
    QGridLayout,
    QHeaderView,
    QSizePolicy,
    QTableWidget,
    QWidget,
)


TABLE_CHECKBOX_COLUMN_WIDTH = 52
TABLE_CHECKBOX_SIZE = 28
TABLE_CHECKBOX_ROW_MIN_HEIGHT = 40
TABLE_CHECKBOX_INDICATOR_SIZE = 18


class CenteredTableCheckBox(QCheckBox):
    """QCheckBox variant with a manually centered indicator.

    Qt stylesheets lay out ``QCheckBox::indicator`` differently across
    platforms. The table selection column needs pixel-stable centering, so this
    widget keeps the QCheckBox API but paints the indicator directly.
    """

    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.setAttribute(Qt.WA_Hover, True)

    def sizeHint(self):
        return QSize(TABLE_CHECKBOX_SIZE, TABLE_CHECKBOX_SIZE)

    def paintEvent(self, event):
        indicator_rect = QRect(
            QPoint(0, 0),
            QSize(TABLE_CHECKBOX_INDICATOR_SIZE, TABLE_CHECKBOX_INDICATOR_SIZE),
        )
        indicator_rect.moveCenter(self.rect().center())

        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)

        palette = self.palette()
        checked = self.isChecked()
        partially_checked = self.checkState() == Qt.CheckState.PartiallyChecked
        enabled = self.isEnabled()
        hovered = self.underMouse()

        if checked or partially_checked:
            fill_color = palette.highlight().color()
            border_color = palette.highlight().color()
        else:
            fill_color = palette.base().color() if enabled else palette.button().color()
            border_color = palette.highlight().color() if enabled and hovered else palette.mid().color()

        rect = QRectF(indicator_rect).adjusted(0.5, 0.5, -0.5, -0.5)
        painter.setPen(QPen(border_color, 1))
        painter.setBrush(fill_color)
        painter.drawRoundedRect(rect, 4, 4)

        if checked:
            pen = QPen(palette.highlightedText().color(), 2)
            pen.setCapStyle(Qt.RoundCap)
            pen.setJoinStyle(Qt.RoundJoin)
            painter.setPen(pen)
            painter.drawLine(
                QPointF(rect.left() + 4, rect.center().y() + 1),
                QPointF(rect.left() + 7, rect.bottom() - 4),
            )
            painter.drawLine(
                QPointF(rect.left() + 7, rect.bottom() - 4),
                QPointF(rect.right() - 4, rect.top() + 5),
            )
        elif partially_checked:
            pen = QPen(palette.highlightedText().color(), 2)
            pen.setCapStyle(Qt.RoundCap)
            painter.setPen(pen)
            painter.drawLine(
                QPointF(rect.left() + 4, rect.center().y()),
                QPointF(rect.right() - 4, rect.center().y()),
            )

    def hitButton(self, pos):
        return self.rect().contains(pos)

    def enterEvent(self, event):
        super().enterEvent(event)
        self.update()

    def leaveEvent(self, event):
        super().leaveEvent(event)
        self.update()


class TableCheckBoxCell(QWidget):
    def __init__(self, checkbox: QCheckBox, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.checkbox = checkbox
        self.setObjectName("tableCheckBoxCell")
        self.setMinimumHeight(TABLE_CHECKBOX_ROW_MIN_HEIGHT)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        layout = QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.checkbox, 0, 0, Qt.AlignCenter)


def configure_table_checkbox(checkbox: QCheckBox) -> QCheckBox:
    checkbox.setObjectName("tableCheckBox")
    checkbox.setText("")
    checkbox.setFixedSize(TABLE_CHECKBOX_SIZE, TABLE_CHECKBOX_SIZE)
    checkbox.setFocusPolicy(Qt.NoFocus)
    checkbox.setCursor(Qt.PointingHandCursor)
    return checkbox


def create_table_checkbox(parent: Optional[QWidget] = None) -> QCheckBox:
    return configure_table_checkbox(CenteredTableCheckBox(parent))


def make_table_checkbox_cell(checkbox: QCheckBox) -> QWidget:
    return TableCheckBoxCell(checkbox)


def configure_table_checkbox_column(table: QTableWidget) -> None:
    header = table.horizontalHeader()
    header.setSectionResizeMode(0, QHeaderView.Fixed)
    header.setMinimumHeight(max(header.minimumHeight(), TABLE_CHECKBOX_SIZE + 8))
    table.setColumnWidth(0, TABLE_CHECKBOX_COLUMN_WIDTH)
    vertical_header = table.verticalHeader()
    vertical_header.setSectionResizeMode(QHeaderView.Fixed)
    vertical_header.setDefaultSectionSize(
        max(vertical_header.defaultSectionSize(), TABLE_CHECKBOX_ROW_MIN_HEIGHT)
    )


def position_header_checkbox(table: QTableWidget, checkbox: QCheckBox) -> None:
    header = table.horizontalHeader()
    if not header:
        return
    x = header.sectionPosition(0)
    w = header.sectionSize(0)
    h = header.height()
    size = TABLE_CHECKBOX_SIZE
    checkbox.setGeometry(
        x + max(0, (w - size) // 2),
        max(0, (h - size) // 2),
        size,
        size,
    )
