"""自定义异常与统一错误响应。"""

from __future__ import annotations

from typing import Any, List, Optional


class APIError(Exception):
    def __init__(
        self,
        status_code: int,
        code: str,
        message: str,
        details: Any = None,
    ) -> None:
        super().__init__(message)
        self.status_code = status_code
        self.code = code
        self.message = message
        self.details = details or {}


class ValidationError(APIError):
    def __init__(self, errors: List[dict], message: str = "请求载荷校验失败") -> None:
        super().__init__(422, "validation_failed", message, {"errors": errors})


class PayloadConflictError(APIError):
    def __init__(self, code: str, message: str, details: Optional[dict] = None) -> None:
        super().__init__(409, code, message, details or {})


class NotFound(APIError):
    def __init__(self, code: str, message: str, details: Optional[dict] = None) -> None:
        super().__init__(404, code, message, details or {})
