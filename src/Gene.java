import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class Gene {
    String id;
    String symbol;
    char strand;
    String chromosome;
    HashMap<String, Transcript> transcriptHashMap = new HashMap<>();
    Set<String> proteinHashSet = new HashSet<>();
    HashMap<Intron, Set<Transcript>> uniqueIntronMap = new HashMap<>();
    HashMap<CDS, Set<Transcript>> uniqueCDSMap = new HashMap<>();

    public Gene(String id, String symbol, char strand, String chromosome) {
        this.id = id;
        this.symbol = symbol;
        this.strand = strand;
        this.chromosome = chromosome;
    }
    public void addUniqueIntron(Intron intron, Transcript transcript){
        for (Intron key : uniqueIntronMap.keySet()){
            if (key.equals(intron)){
                uniqueIntronMap.get(key).add(transcript);
                return;
            }
        }
        Set<Transcript> transcriptSet= new HashSet<>();
        transcriptSet.add(transcript);
        uniqueIntronMap.put(intron, transcriptSet);
    }
    public void addUniqueCDS(CDS cds, Transcript transcript){
        for (CDS key : uniqueCDSMap.keySet()){
            if (key.equals(cds)){
                uniqueCDSMap.get(key).add(transcript);
                return;
            }
        }
        Set<Transcript> transcriptSet= new HashSet<>();
        transcriptSet.add(transcript);
        uniqueCDSMap.put(cds, transcriptSet);
    }
}
