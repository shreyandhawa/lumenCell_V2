"""Error handling: every failure leaves the server as an ``ErrorBody``."""

from http import HTTPStatus
from typing import Any

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from starlette.exceptions import HTTPException as StarletteHTTPException

from lumencell.api.models import ErrorBody, ErrorInfo


class ApiError(Exception):
    """Raise from a route to send a structured error with a stable code."""

    def __init__(self, status: int, code: str, message: str, detail: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.status = status
        self.code = code
        self.message = message
        self.detail = detail


def root_required() -> ApiError:
    return ApiError(
        HTTPStatus.CONFLICT,
        "root_required",
        "Pick where the trajectory starts before running a method. "
        "lumenCell never guesses the root, because a wrong guess silently reverses pseudotime.",
    )


def _body(code: str, message: str, detail: dict[str, Any] | None = None) -> dict[str, Any]:
    return ErrorBody(error=ErrorInfo(code=code, message=message, detail=detail)).model_dump()


def install_error_handlers(app: FastAPI) -> None:
    @app.exception_handler(ApiError)
    async def _api_error(_request: Request, exc: ApiError) -> JSONResponse:
        return JSONResponse(_body(exc.code, exc.message, exc.detail), status_code=exc.status)

    @app.exception_handler(RequestValidationError)
    async def _validation(_request: Request, exc: RequestValidationError) -> JSONResponse:
        errors = [
            {"loc": [str(p) for p in e.get("loc", ())], "msg": str(e.get("msg", ""))} for e in exc.errors()
        ]
        return JSONResponse(
            _body("invalid_request", "The request is not valid.", {"errors": errors}),
            status_code=HTTPStatus.UNPROCESSABLE_ENTITY,
        )

    @app.exception_handler(StarletteHTTPException)
    async def _http(_request: Request, exc: StarletteHTTPException) -> JSONResponse:
        code = {404: "not_found", 405: "method_not_allowed"}.get(exc.status_code, "http_error")
        return JSONResponse(_body(code, str(exc.detail)), status_code=exc.status_code)

    @app.exception_handler(Exception)
    async def _unexpected(_request: Request, _exc: Exception) -> JSONResponse:
        # Anything a route did not anticipate. The message stays generic (the
        # exception text may hold paths); the server re-raises after this
        # response is sent, so the traceback is in the terminal running lumencell.
        return JSONResponse(
            _body(
                "internal_error",
                "Something went wrong on the server. The terminal running lumencell has the details.",
            ),
            status_code=HTTPStatus.INTERNAL_SERVER_ERROR,
        )
