import java.util.Objects;

public class Intron extends CDS {

    public Intron(Transcript source, int start, int end, char strand) {
        super(source, "", start, end, strand);
    }

    public boolean contains(Intron intron) {
        return (this.start <= intron.start && this.end >= intron.end && this.strand == intron.strand);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true; // same object
        if (!(o instanceof Intron)) return false; // wrong type
        Intron i = (Intron) o;
        return start == i.start &&
                end == i.end &&
                strand == i.strand;
    }
}