from typing import TypeVar, Generic, List, Optional, Union
import  enum

T = TypeVar('T')


class Annotation(Generic[T]):
    """
    Generic annotation class to store origin and data of any type.

    Args:
        origin: String identifier for the annotation source
        data: The annotation data of type T
    """

    def __init__(self, origin: str, data: T) -> None:
        self._origin = origin
        self._data = data

    @property
    def origin(self) -> str:
        """Get the annotation origin."""
        return self._origin

    @origin.setter
    def origin(self, value: str) -> None:
        """Set the annotation origin."""
        if not isinstance(value, str):
            raise TypeError("Origin must be a string")
        self._origin = value

    @property
    def data(self) -> T:
        """Get the annotation data."""
        return self._data

    @data.setter
    def data(self, value: T) -> None:
        """Set the annotation data."""
        self._data = value


class Annotation_Float (Annotation[float]):
    """Annotation specifically for float data."""

    def __init__(self, origin: str, data: float):
        if not isinstance(data, (int, float)):
            raise TypeError("Data must be a number")
        super().__init__(origin, float(data))

    def round_data(self, decimals: int = 2) -> float:
        """Round the float data to specified decimal places."""
        return round(self.data, decimals)


class Annotation_Enum(Annotation[enum.Enum]):
    """Annotation for enumeration data."""

    def __init__(self, origin: str, data: enum.Enum):
        if not isinstance(data, enum.Enum):
            raise TypeError("Data must be an Enum")
        super().__init__(origin, data)

class Annotation_Str(Annotation[Union[str, List[str]]]):
    """Annotation for string or list of strings data."""

    def __init__(self, origin: str, data: Union[str, List[str]]):
        if not isinstance(data, (str, list)):
            raise TypeError("Data must be a string or list of strings")
        if isinstance(data, list) and not all(isinstance(item, str) for item in data):
            raise TypeError("All items in list must be strings")
        super().__init__(origin, data)

    def get_as_string(self, separator: str = ", ") -> str:
        """Convert data to string, joining list items if necessary."""
        if isinstance(self.data, list):
            return separator.join(self.data)
        return self.data

    def get_length(self) -> int:
        """Get length of data (string length or list length)."""
        return len(self.data)


class TranscriptAnnotations:
    """Standalone class to manage transcript-specific annotation values"""

    def __init__(self, transcript_id: str = "", impact: str="", revel: float = 0.0, gerp: float = 0.0, spliceai: float = 0.0,
                 polyphen: float = 0.0, loftee: str = "",  is_mane: bool = False,
                 cdna_notation: str = "", protein_notation: str = "", consequences: str = "") -> None:
        self.transcript_id = transcript_id
        self.impact = impact
        self._revel = revel
        self._gerp = gerp
        self._spliceai = spliceai
        self._polyphen = polyphen
        self._loftee = loftee
        self._is_mane = is_mane
        self._cdna_notation = cdna_notation
        self._protein_notation = protein_notation
        self._consequences = consequences

    def add_values(self, transcript_id: str = None, impact: str = None,  revel: float = None, gerp: float = None, spliceai: float = None,
                   polyphen: float = None, loftee: str = None, is_mane: bool = None,
                   cdna_notation: str = None, protein_notation: str = None, consequences: str = None) -> None:
        """Add or update annotation values"""
        if transcript_id is not None:
            self._transcript_id = transcript_id
        if impact is not None:
            self._impact = impact
        if revel is not None:
            self._revel = revel
        if gerp is not None:
            self._gerp = gerp
        if spliceai is not None:
            self._spliceai = spliceai
        if polyphen is not None:
            self._polyphen = polyphen
        if loftee is not None:
            self._loftee = loftee
        if is_mane is not None:
            self._is_mane = is_mane
        if cdna_notation is not None:
            self._cdna_notation = cdna_notation
        if protein_notation is not None:
            self._protein_notation = protein_notation
        if consequences is not None:
            self._consequences = consequences



class GeneAnnotations:
    """Class to manage a collection of gene annotations"""

    def __init__(self, variant_annotations: List[Annotation] = None,  transcript_annotations: List[TranscriptAnnotations]= None) -> None:
        """Initialize with a list of annotations"""
        self._variant_annotations = variant_annotations or []
        self._transcript_annotations = transcript_annotations or []

    def add_annotation(self, annotation: Annotation) -> None:
        """Add a single variant specific annotation to the collection"""
        if annotation is not None:
            self._variant_annotations.append(annotation)

    def add_transcript_annotation(self, transcript_annotation: TranscriptAnnotations) -> None:
        """Add a transcript annotation to the collection"""
        if transcript_annotation is not None:
            self._transcript_annotations.append(transcript_annotation)

    def get_all_annotations(self) -> List[Annotation]:
        """Get all variant specific annotations"""
        return self._variant_annotations.copy()

    def get_transcript_annotations(self) -> List[TranscriptAnnotations]:
        """Get all transcript annotations"""
        return self._transcript_annotations.copy()

    def get_annotation_by_origin(self, origin: str) -> Optional[Annotation]:
        """Get first annotation with specified origin, returns None if not found"""
        for ann in self._variant_annotations:
            if ann.origin == origin:
                return ann
        return None

    def has_origin(self, origin: str) -> bool:
        """Check if any annotation has the specified origin"""
        return any(ann.origin == origin for ann in self._variant_annotations)


class ImpactLevel(enum.Enum):
    """VEP impact levels"""
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MODIFIER = "MODIFIER"

class LofteeLevel(enum.Enum):
    """LoF (Loss of Function) confidence levels"""
    HC = "HC"  # High Confidence
    LC = "LC"  # Low Confidence

