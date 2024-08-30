from sqlalchemy import Column, Integer, String, DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.sql import func

Base = declarative_base()

class Molecule(Base):
    __tablename__ = "molecule"
    id = Column(Integer, primary_key=True, autoincrement=True)
    identifier = Column(String)
    smiles = Column(String)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())

    # def __repr__(self):
    #     return f"<Molecule(id={self.id}, identifier={self.identifier}, smiles={self.smiles})>"
    
    # def __str__(self):
    #     return f"<Molecule(id={self.id}, identifier={self.identifier}, smiles={self.smiles})>"