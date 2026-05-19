"""Shared layout helpers for analysis input pages."""
from typing import Tuple, Type

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QFrame,
    QGroupBox,
    QLayout,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QTextEdit,
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


def configure_pairwise_limit_controls(
    candidate_limit_input: QSpinBox,
    pairwise_all_input: QCheckBox,
    default_value: int = 3,
) -> None:
    """Keep Pairwise Top-N and full-candidate mode semantically aligned."""
    normal_min = candidate_limit_input.minimum()
    normal_max = candidate_limit_input.maximum()
    default_value = min(max(default_value, normal_min), normal_max)
    candidate_limit_input.setProperty("lastPairwiseLimit", candidate_limit_input.value() or default_value)
    candidate_limit_input.setToolTip("限制每个查询基因组参与 pairwise 的候选数量")
    pairwise_all_input.setToolTip("勾选后使用全部候选，忽略 Pairwise Top-N")

    def sync_candidate_limit(checked: bool) -> None:
        if checked:
            current_value = candidate_limit_input.value()
            if normal_min <= current_value <= normal_max:
                candidate_limit_input.setProperty("lastPairwiseLimit", current_value)
            candidate_limit_input.setRange(0, 0)
            candidate_limit_input.setSpecialValueText(" ")
            candidate_limit_input.setValue(0)
            candidate_limit_input.setEnabled(False)
            candidate_limit_input.setToolTip("全量候选模式下不使用 Pairwise Top-N 限制")
            return

        saved_value = candidate_limit_input.property("lastPairwiseLimit") or default_value
        try:
            restored_value = int(saved_value)
        except (TypeError, ValueError):
            restored_value = default_value
        restored_value = min(max(restored_value, normal_min), normal_max)
        candidate_limit_input.setEnabled(True)
        candidate_limit_input.setSpecialValueText("")
        candidate_limit_input.setRange(normal_min, normal_max)
        candidate_limit_input.setValue(restored_value)
        candidate_limit_input.setToolTip("限制每个查询基因组参与 pairwise 的候选数量")

    pairwise_all_input.toggled.connect(sync_candidate_limit)
    sync_candidate_limit(pairwise_all_input.isChecked())


def text_edit_height_for_lines(text_edit: QTextEdit, line_count: int, min_lines: int = 2) -> int:
    """Calculate a compact QTextEdit height from its visible line count."""
    return text_edit.fontMetrics().lineSpacing() * max(min_lines, line_count) + 30


def configure_auto_growing_text_edit(text_edit: QTextEdit, min_lines: int = 2):
    """Grow a QTextEdit with its content while leaving page scroll to the parent."""
    text_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    text_edit.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

    def resize_to_content() -> None:
        line_count = text_edit.toPlainText().count("\n") + 1
        text_edit.setFixedHeight(text_edit_height_for_lines(text_edit, line_count, min_lines))
        text_edit.updateGeometry()

    text_edit.textChanged.connect(resize_to_content)
    resize_to_content()
    return resize_to_content
