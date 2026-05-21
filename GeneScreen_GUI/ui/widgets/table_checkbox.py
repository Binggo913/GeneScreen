"""Shared table checkbox helpers for selectable table columns."""

from typing import Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox, QHeaderView, QSizePolicy, QTableWidget, QWidget


TABLE_CHECKBOX_COLUMN_WIDTH = 52
TABLE_CHECKBOX_SIZE = 28
TABLE_CHECKBOX_ROW_MIN_HEIGHT = 40


class TableCheckBoxCell(QWidget):
    def __init__(self, checkbox: QCheckBox, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.checkbox = checkbox
        self.checkbox.setParent(self)
        self.setObjectName("tableCheckBoxCell")
        self.setMinimumHeight(TABLE_CHECKBOX_ROW_MIN_HEIGHT)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._center_checkbox()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._center_checkbox()

    def _center_checkbox(self) -> None:
        self.checkbox.move(
            max(0, (self.width() - self.checkbox.width()) // 2),
            max(0, (self.height() - self.checkbox.height()) // 2),
        )


def configure_table_checkbox(checkbox: QCheckBox) -> QCheckBox:
    checkbox.setObjectName("tableCheckBox")
    checkbox.setText("")
    checkbox.setFixedSize(TABLE_CHECKBOX_SIZE, TABLE_CHECKBOX_SIZE)
    checkbox.setFocusPolicy(Qt.NoFocus)
    checkbox.setCursor(Qt.PointingHandCursor)
    return checkbox


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
