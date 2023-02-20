from __future__ import annotations

from typing import Final

import sqlalchemy
from sqlalchemy import Column

from libysjs.ds.ysjs_job import YSJSJob
from ysjsd.orm import SQLAlchemyDeclarativeBase


class YSJYSJSJobTable(SQLAlchemyDeclarativeBase):
    __tablename__: Final[str] = "ysjs_job"
    job_id = Column(sqlalchemy.Integer, primary_key=True, nullable=False)
    submission_id = Column(sqlalchemy.String(32), nullable=False)
    status = Column(sqlalchemy.String(8), nullable=False)
    retv = Column(sqlalchemy.Integer, nullable=True)
    start_time = Column(sqlalchemy.Float, nullable=True)
    terminate_time = Column(sqlalchemy.Float, nullable=True)
    pid = Column(sqlalchemy.Integer, nullable=True)

    @classmethod
    def from_job(cls, job: YSJSJob):
        submission_id = job.submission.submission_id
        job_dict = dict(job.to_dict())
        _ = job_dict.pop("submission")
        return cls(
            submission_id=submission_id,
            **job_dict
        )
