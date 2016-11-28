#pragma once

// ----------------------------------------------------------------------

enum class SeqdbJsonKey : char
{
    Comment='?', Comment_='_',
    Name='N', Dates='d', Continent='C', Country='c', Lineage='l', VirusType='v',
    SequenceSet='s',
    AminoAcids='a', Nucleotides='n', Clades='c', Gene='g', HiNames='h', LabIds='l',
    Passages='p', Reassortant='r', AminoAcidShift='s', NucleotideShift='t',

    Unknown
};

// ----------------------------------------------------------------------
/// Local Variables:
/// eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
/// End:
