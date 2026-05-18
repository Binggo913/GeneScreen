"""Shared layout helpers for analysis input pages."""
from typing import Tuple, Type

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame,
    QGroupBox,
    QLayout,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


CARD_MARGINS: Tuple[int, int, int, int] = (14, 24, 14, 14)


def create_scroll_content(
    parent: QWidget,
    margins: Tuple[int, int, int, int] = (30, 30, 30, 30),
    spacing: int = 20,
) -> QVBoxLayout:
    root_layout = QVBoxLayout(parent)
    root_layout.setContentsMargins(0, 0, 0, 0)
    root_layout.setSpacing(0)

    scroll_area = QScrollArea(parent)
    scroll_area.setWidgetResizable(True)
    scroll_area.setFrameShape(QFrame.NoFrame)
    scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

    content = QWidget()
    content_layout = QVBoxLayout(content)
    content_layout.setContentsMargins(*margins)
    content_layout.setSpacing(spacing)

    scroll_area.setWidget(content)
    root_layout.addWidget(scroll_area)
    return content_layout


def create_card(
    title: str,
    layout_cls: Type[QLayout] = QVBoxLayout,
    margins: Tuple[int, int, int, int] = CARD_MARGINS,
    spacing: int = 12,
) -> Tuple[QGroupBox, QLayout]:
    card = QGroupBox(title)
    card.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)

    layout = layout_cls(card)
    layout.setContentsMargins(*margins)
    layout.setSpacing(spacing)
    return card, layout
