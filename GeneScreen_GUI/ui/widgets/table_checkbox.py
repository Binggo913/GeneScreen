"""Shared table checkbox helpers for selectable table columns."""

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox, QGridLayout, QHeaderView, QTableWidget, QWidget


TABLE_CHECKBOX_COLUMN_WIDTH = 52
TABLE_CHECKBOX_SIZE = 28
TABLE_CHECKBOX_ROW_MIN_HEIGHT = 40


def configure_table_checkbox(checkbox: QCheckBox) -> QCheckBox:
    checkbox.setObjectName("tableCheckBox")
    checkbox.setText("")
    checkbox.setFixedSize(TABLE_CHECKBOX_SIZE, TABLE_CHECKBOX_SIZE)
    checkbox.setFocusPolicy(Qt.NoFocus)
    checkbox.setCursor(Qt.PointingHandCursor)
    return checkbox


def make_table_checkbox_cell(checkbox: QCheckBox) -> QWidget:
    cell = QWidget()
    cell.setObjectName("tableCheckBoxCell")
    cell.setMinimumHeight(TABLE_CHECKBOX_ROW_MIN_HEIGHT)
    layout = QGridLayout(cell)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(0)
    layout.addWidget(checkbox, 0, 0, Qt.AlignCenter)
    return cell


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
