from __future__ import annotations

from typing import Final

import sqlalchemy
from sqlalchemy import Column

from ysjsd.orm import SQLAlchemyDeclarativeBase


class ServerSideYSJSDConfigTable(SQLAlchemyDeclarativeBase):
    __tablename__: Final[str] = "ysjsd_config"
    name = Column(sqlalchemy.String(32), primary_key=True, nullable=False)
    description = Column(sqlalchemy.String(1024), nullable=False)
    ysjs_port = Column(sqlalchemy.String(8), nullable=False)
    var_directory_path = Column(sqlalchemy.String(256), nullable=False)
    config_file_path = Column(sqlalchemy.String(256), nullable=False)
    total_cpu = Column(sqlalchemy.Float, nullable=False)
    total_mem = Column(sqlalchemy.Float, nullable=False)
    schedule_method = Column(sqlalchemy.String(8), nullable=False)
    max_concurrent_jobs = Column(sqlalchemy.Integer, nullable=False)
    kill_timeout = Column(sqlalchemy.Float, nullable=False)
