from __future__ import annotations

from typing import Final

import sqlalchemy
from sqlalchemy import Column

from ysjsd.orm import SQLAlchemyDeclarativeBase


class YSJSSubmissionTable(SQLAlchemyDeclarativeBase):
    __tablename__: Final[str] = "ysjs_submission"
    submission_id = Column(sqlalchemy.String(32), primary_key=True, nullable=False)
    submission_name = Column(sqlalchemy.String(32), nullable=False)
    submission_description = Column(sqlalchemy.String(4096), nullable=False)
    cpu = Column(sqlalchemy.Float, nullable=False)
    mem = Column(sqlalchemy.Float, nullable=False)
    submission_time = Column(sqlalchemy.Float, nullable=False)
    cwd = Column(sqlalchemy.String(256), nullable=False)
    tags = Column(sqlalchemy.JSON, nullable=False)
    env = Column(sqlalchemy.JSON, nullable=False)
    stdin = Column(sqlalchemy.String(256), nullable=True)
    stdout = Column(sqlalchemy.String(256), nullable=True)
    stderr = Column(sqlalchemy.String(256), nullable=True)
    script_path = Column(sqlalchemy.String(256), nullable=False)
    shell_path = Column(sqlalchemy.String(256), nullable=False)
    depends = Column(sqlalchemy.JSON, nullable=False)
