from __future__ import annotations

from typing import Final

import sqlalchemy
from sqlalchemy import Column
from sqlalchemy.orm import declarative_base

SQLAlchemyDeclarativeBase = declarative_base()


class ServerSideYSJSDConfigTable(SQLAlchemyDeclarativeBase):
    __tablename__: Final[str] = "ysjsd_config"
    name = Column(sqlalchemy.String(32), primary_key=True)
    description = Column(sqlalchemy.String(1024))
    ysjs_port = Column(sqlalchemy.String(8))
    var_directory = Column(sqlalchemy.String(256))
    config_file_path = Column(sqlalchemy.String(256))
    total_cpu = Column(sqlalchemy.Float)
    total_mem = Column(sqlalchemy.Float)
    schedule_method = Column(sqlalchemy.String(8))
    max_concurrent_jobs = Column(sqlalchemy.Integer)
    kill_timeout = Column(sqlalchemy.Float)



class YSJSDVersionTable(SQLAlchemyDeclarativeBase):
    __tablename__: Final[str] = "ysjsd_versions"
    name = Column(sqlalchemy.String(32), primary_key=True)
    version = Column(sqlalchemy.String(32))


class SubmissionTable(SQLAlchemyDeclarativeBase):
    __tablename__: Final[str] = "ysjsd_submissions"
    submission_id = Column(sqlalchemy.String(32), primary_key=True)
    submission_name = Column(sqlalchemy.String(32))
    submission_description = Column(sqlalchemy.String(4096))
    cpu = Column(sqlalchemy.Float)
    mem = Column(sqlalchemy.Float)
    submission_time = Column(sqlalchemy.Float)
    cwd = Column(sqlalchemy.String(256))
    tags = Column(sqlalchemy.JSON)
    env = Column(sqlalchemy.JSON)
    stdin = Column(sqlalchemy.String(256))
    stdout = Column(sqlalchemy.String(256))
    stderr = Column(sqlalchemy.String(256))
    script_path = Column(sqlalchemy.String(256))
    shell_path = Column(sqlalchemy.String(256))
