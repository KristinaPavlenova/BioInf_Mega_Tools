import pytest
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from bioinf_mega_tools import filter_fastq, parse_range


@pytest.fixture
def example_fastq(tmp_path):
    records = [
        SeqRecord(Seq("ATGC"), id="read1", letter_annotations={"phred_quality": [40, 40, 40, 40]}),
        SeqRecord(Seq("GCGC"), id="read2", letter_annotations={"phred_quality": [10, 10, 10, 10]}),
        SeqRecord(Seq("AAAAA"), id="read3", letter_annotations={"phred_quality": [30, 30, 30, 30, 30]}),
    ]
    fastq_path = tmp_path / "test.fastq"
    with open(fastq_path, "w") as handle:
        SeqIO.write(records, handle, "fastq")
    return fastq_path


class TestFilterFastq:
    def test_quality_filtering(self, example_fastq, tmp_path):
        output_path = tmp_path / "out.fastq"
        filter_fastq(str(example_fastq), str(output_path), quality_threshold=20)
        records = list(SeqIO.parse(output_path, "fastq"))
        assert all(sum(record.letter_annotations["phred_quality"]) / len(record.seq) >= 20 for record in records)


    def test_gc_bounds_filtering(self, example_fastq, tmp_path):
        output_path = tmp_path / "out.fastq"
        filter_fastq(str(example_fastq), str(output_path), gc_bounds=(60, 100))
        records = list(SeqIO.parse(output_path, "fastq"))
        assert all(60 <= (record.seq.count("G") + record.seq.count("C")) / len(record.seq) * 100 <= 100 for record in records)


    def test_length_bounds_filtering(self, example_fastq, tmp_path):
        output_path = tmp_path / "out.fastq"
        filter_fastq(str(example_fastq), str(output_path), length_bounds=(5, 10))
        records = list(SeqIO.parse(output_path, "fastq"))
        assert all(5 <= len(record.seq) <= 10 for record in records)


    def test_output_file_created(self, example_fastq, tmp_path):
        output_path = tmp_path / "out.fastq"
        filter_fastq(str(example_fastq), str(output_path))
        assert output_path.exists()


class TestParseRange:
    def test_parse_range_tuple(self):
        result = parse_range(None, None, "20,80")
        assert result == (20.0, 80.0)


    def test_parse_range_single(self):
        result = parse_range(None, None, "75")
        assert result == 75.0


class TestErrorHandling:
    def test_nonexistent_file(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            filter_fastq("non_existent_file.fastq", str(tmp_path / "out.fastq"))


class TestEmptyFile:
    def test_empty_fastq(self, tmp_path):
        empty_file = tmp_path / "empty.fastq"
        empty_file.write_text("")
        output_path = tmp_path / "out.fastq"
        filter_fastq(str(empty_file), str(output_path))
        assert os.path.exists(output_path)
        assert list(SeqIO.parse(output_path, "fastq")) == []
