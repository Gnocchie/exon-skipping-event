import java.util.Objects;

public class CDS {
    Transcript source;
    String prot_id;
    int start;
    int end;
    char strand;


    public CDS(Transcript source, String prot_id, int start, int end, char strand) {
        this.source = source;
        this.prot_id = prot_id;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }
    public boolean equals(Object o){
        if (this == o) return true; // same object
        if (!(o instanceof CDS)) return false; // wrong type
        CDS cds = (CDS) o;
        return start == cds.start &&
                end == cds.end &&
                strand == cds.strand;    }

    @Override
    public int hashCode() {
        return Objects.hash(start, end, strand);
    }
    public boolean contains(CDS cds){
        return (this.start <= cds.start && this.end >= cds.end && this.strand == cds.strand);
    }
}
