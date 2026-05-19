"""Global background analysis task queue for the GUI."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from typing import Callable, Deque, List, Optional

from PySide6.QtCore import QObject, QThread, Signal

from .config import get_analysis_thread_count
from .database import get_database


@dataclass
class AnalysisTask:
    thread: QThread
    history_ids: List[int]
    label: str = ""
    on_started: Optional[Callable[[], None]] = None


class AnalysisTaskManager(QObject):
    """Run submitted analysis QThreads with a configurable concurrency limit."""

    task_counts_changed = Signal(int, int)  # running, pending

    def __init__(self, parent=None):
        super().__init__(parent)
        self._max_running = get_analysis_thread_count()
        self._pending: Deque[AnalysisTask] = deque()
        self._running: List[AnalysisTask] = []

    def set_max_running(self, count: int) -> None:
        self._max_running = max(1, min(int(count), 32))
        self._start_available()
        self._emit_counts()

    def submit(self, task: AnalysisTask) -> None:
        self._pending.append(task)
        self._start_available()
        self._emit_counts()

    def running_count(self) -> int:
        return len(self._running)

    def pending_count(self) -> int:
        return len(self._pending)

    def has_tasks(self) -> bool:
        return bool(self._running or self._pending)

    def mark_unfinished_failed(self) -> None:
        for task in list(self._running) + list(self._pending):
            for history_id in task.history_ids:
                get_database().update_history(history_id, status="failed")

    def _start_available(self) -> None:
        while self._pending and len(self._running) < self._max_running:
            task = self._pending.popleft()
            self._running.append(task)
            for history_id in task.history_ids:
                get_database().update_history(history_id, status="running")
            if task.on_started:
                task.on_started()
            task.thread.finished.connect(lambda *_, t=task: self._on_task_finished(t))
            task.thread.start()

    def _on_task_finished(self, task: AnalysisTask) -> None:
        if task in self._running:
            self._running.remove(task)
        task.thread.deleteLater()
        self._start_available()
        self._emit_counts()

    def _emit_counts(self) -> None:
        self.task_counts_changed.emit(self.running_count(), self.pending_count())


_task_manager: Optional[AnalysisTaskManager] = None


def get_analysis_task_manager() -> AnalysisTaskManager:
    global _task_manager
    if _task_manager is None:
        _task_manager = AnalysisTaskManager()
    return _task_manager
