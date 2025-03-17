from pydantic import BaseModel  # pydantic 2.0+


class Msg(BaseModel):
    """Schema for generic message responses in API endpoints."""
    msg: str