import java.util.*;

public class Transcript {
    String id;
    Gene source;
    List<CDS> cdsList = new ArrayList<>();
    List<Intron> intronList = new ArrayList<>();
    Set<String> proteinList = new HashSet<>();

    public Transcript(String id, Gene source) {
        this.id = id;
        this.source = source;
    }
    public boolean equals(Object o){
        if (this == o) return true; // same object
        if (!(o instanceof Transcript)) return false; // wrong type
        Transcript tr = (Transcript) o;
        return Objects.equals(tr.id, this.id);
    }
    public int hashCode() {
        return Objects.hash(id);
    }
    public CDS getCDSBefore(int pos) {
        CDS prev = null;
        for (CDS cds : cdsList) {
            if (cds.end < pos) {
                if (prev == null || cds.end > prev.end) {
                    prev = cds;
                }
            }
        }
        return prev;
    }

    public CDS getCDSAfter(int pos) {
        CDS next = null;
        for (CDS cds : cdsList) {
            if (cds.start > pos) {
                if (next == null || cds.start < next.start) {
                    next = cds;
                }
            }
        }
        return next;
    }

    public List<CDS> getCDSInRange(int start, int end) {
        List<CDS> inside = new ArrayList<>();
        for (CDS cds : cdsList) {
            if (cds.start >= start && cds.end <= end) {
                inside.add(cds);
            }
        }
        return inside;
    }

}
