from __future__ import annotations

import sqlalchemy
from sqlalchemy import Column

from labw_utils.typing_importer import Final
from ysjsd.orm import SQLAlchemyDeclarativeBase


class YSJSDVersionTable(SQLAlchemyDeclarativeBase):
    __tablename__: Final[str] = "ysjsd_versions"
    name = Column(sqlalchemy.String(32), primary_key=True, nullable=False)
    version = Column(sqlalchemy.String(32), nullable=False)
