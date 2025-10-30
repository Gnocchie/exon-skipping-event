import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class ESSE_Finder {
    ArrayList<ESSE> eventList = new ArrayList<>();
    HashMap<String, Gene> geneHashMap = new HashMap<>();
    public void findEvents(String file) throws IOException {
        ReaderGTF.read(file, geneHashMap);
        calculateEvents();
        printEvents();
    }
    public void findEvents(String gtf, String outPath) throws IOException {
        ReaderGTF.read(gtf, geneHashMap);
        calculateEvents();
        writeEvents(outPath);
    }
    public void calculateEvents() {
        sortCDS();
        calculateIntrons();
        printGeneOverview();
        compareIntrons();
    }
    public void sortCDS(){
        for (Gene gene : geneHashMap.values()) {
            for (Transcript transcript : gene.transcriptHashMap.values()) {
                transcript.cdsList.sort(Comparator.comparingInt(c -> c.start));
            }
        }
    }
    public void calculateIntrons(){
        for (Gene gene : geneHashMap.values()) {
            for (Transcript transcript : gene.transcriptHashMap.values()) {
                List<CDS> cdsList = transcript.cdsList;
                // iterate through cds for introns
                CDS pointer = null;
                for (CDS cds : cdsList){
                    if (pointer==null) {
                        pointer = cds;
                        continue;
                    }
                    if (pointer.end < cds.start){
                        Intron intron = new Intron(transcript, pointer.end + 1, cds.start, cds.strand);
                        transcript.intronList.add(intron);
                        gene.addUniqueIntron(intron, transcript);
                    }
                    pointer = cds;
                }
            }
        }
    }
    public void compareIntrons(){
        for (Gene g : geneHashMap.values()) {
            for (Map.Entry<Intron, Set<Transcript>> entry : g.uniqueIntronMap.entrySet()){
                Intron svi = entry.getKey();
                Set<Transcript> svTranscripts = entry.getValue();

                ESSE currentEvent = new ESSE(g, svi);

                //obtain all sv proteins
                for (Transcript svT : svTranscripts){
                    currentEvent.SV_prots.addAll(svT.proteinList);
                }

                //obtain WT transcripts / introns / proteins
                for (Transcript wtT : g.transcriptHashMap.values()){
                    int skippedExons = 0;
                    int skippedBases = 0;

                    boolean leadingFound = false;
                    boolean skippedFound = false;
                    boolean trailingFound = false;

                    CDS pointer = null;
                    List<CDS> wt_exons = new ArrayList<>();
                    List<Intron> wt_introns = new ArrayList<>();
                    Set<String> wt_prots = new HashSet<>();

                    if (svTranscripts.contains(wtT)) continue;

                    // go through cds list to check transcript for wildtype
                    for (CDS wtcds : wtT.cdsList){
                        // check cds position
                        if (!leadingFound){
                            if (wtcds.end == svi.start-1) { // leading found
                                leadingFound = true;
                                pointer = wtcds;
                                continue;
                            }
                            if (wtcds.end > svi.start) { // sv_intron passed without leading
                                break;
                            }
                            else continue;
                        }

                        Intron preIntron = new Intron(wtT, pointer.end + 1, wtcds.start, wtcds.strand);
                        if (svi.contains(preIntron)) wt_introns.add(preIntron);

                        // if inside
                        if (svi.contains(wtcds)){
                            skippedFound = true; // there is a match
                            wt_exons.add(wtcds); // add current cds to exon list
                            wt_prots.add(wtcds.prot_id);
                            skippedExons++;
                            skippedBases += wtcds.end - wtcds.start + 1;

                            pointer = wtcds;
                        }
                        else if (wtcds.start == svi.end){ // not inside but trailing
                            trailingFound = true;
                            break;
                        }
                        else break; // not inside and not trailing -> sv intron has been passed
                    } // end of cds loop

                    if (leadingFound && trailingFound && skippedFound) { // validate for border-cds and event

                        currentEvent.addAllWTExon(wt_exons);
                        currentEvent.addAllWTIntron(wt_introns);
                        currentEvent.WT_prots.addAll(wt_prots);
                        currentEvent.setSkippedBases(skippedBases);
                        currentEvent.setSkippedExon(skippedExons);
                    }
                }   // end of transcript loop

                if (!currentEvent.WT_exons.isEmpty()){
                    eventList.add(currentEvent);
                }
            }
        }
    }
    public void printGeneOverview() {
        int totalGeneCount=0;
        int totalTranscriptCount=0;
        int totalCDSCount=0;
        int totalProteinCount = 0;
        for (Map.Entry<String, Gene> geneEntry : this.geneHashMap.entrySet()) {
            Gene gene = geneEntry.getValue();
            totalGeneCount++;
            totalProteinCount += gene.proteinHashSet.size();
            System.out.println("Gene ID: " + gene.id
                    + ", Chrom: " + gene.chromosome
                    + ", Strand: " + gene.strand
                    + ", Transcripts: " + gene.transcriptHashMap.size()
                    + ", Proteins: " + gene.proteinHashSet.size());

            for (Map.Entry<String, Transcript> txEntry : gene.transcriptHashMap.entrySet()) {
                totalTranscriptCount++;
                Transcript tx = txEntry.getValue();
                System.out.print("  Transcript ID: " + tx.id
                        + ", CDS count: " + tx.cdsList.size()
                        + ", Intron count: " + tx.intronList.size()
                        + ", Protein count: " + tx.proteinList.size()
                        + ", CDS positions: ...");
                System.out.println(); // new line after CDS positions
            }
            System.out.println("------------------------------------------------");
        }
        System.out.println("Gene count: " + totalGeneCount);
        System.out.println("Transcript count: " + totalTranscriptCount);
        System.out.println("CDS count: " + totalCDSCount);
        System.out.println("Protein Count: " + totalProteinCount);
    }
    public void printEvents(){
        System.out.println(String.join("\t",
                "id",
                "symbol",
                "chr",
                "strand",
                "nprots",
                "ntrans",
                "SV",
                "WT",
                "WT_prots",
                "SV_prots",
                "min_skipped_exon",
                "max_skipped_exon",
                "min_skipped_bases",
                "max_skipped_bases"));
        for (ESSE event : eventList){

            System.out.println(String.join("\t",
                    event.gene_id,
                    event.gene_symbol,
                    event.chromosome,
                    String.valueOf(event.strand),
                    String.valueOf(event.gene.proteinHashSet.size()),
                    String.valueOf(event.gene.transcriptHashMap.size()),
                    event.getInterval(),
                    event.getWTIntrons(),
                    event.getProtsWT(),
                    event.getProtsSV(),
                    String.valueOf(event.min_skipped_exon),
                    String.valueOf(event.max_skipped_exon),
                    String.valueOf(event.min_skipped_bases),
                    String.valueOf(event.max_skipped_bases)
            ));
        }
    }
    public void writeEvents(String outputPath) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath))) {

            // Write header
            writer.write(String.join("\t",
                    "id",
                    "symbol",
                    "chr",
                    "strand",
                    "nprots",
                    "ntrans",
                    "SV",
                    "WT",
                    "SV_prots",
                    "WT_prots",
                    "min_skipped_exon",
                    "max_skipped_exon",
                    "min_skipped_bases",
                    "max_skipped_bases"));
            writer.newLine();

            // Write each event line
            for (ESSE event : eventList) {
                writer.write(String.join("\t",
                        event.gene_id,
                        event.gene_symbol,
                        event.chromosome,
                        String.valueOf(event.strand),
                        String.valueOf(event.gene.proteinHashSet.size()),
                        String.valueOf(event.gene.transcriptHashMap.size()),
                        event.getInterval(),
                        event.getWTIntrons(),
                        event.getProtsSV(),
                        event.getProtsWT(),
                        String.valueOf(event.min_skipped_exon),
                        String.valueOf(event.max_skipped_exon),
                        String.valueOf(event.min_skipped_bases),
                        String.valueOf(event.max_skipped_bases)
                ));
                writer.newLine();
            }
        } catch (IOException e) {
            System.err.println("Error writing events to file: " + e.getMessage());
        }
    }

}
class ESSE {
    Gene gene;
    String gene_id;
    String gene_symbol;
    String chromosome;
    char strand;
    int svStart;
    int svEnd;
    Set<Intron> WT_introns = new HashSet<>();
    Set<CDS> WT_exons = new HashSet<>();
    Set<String> WT_prots = new HashSet<>();
    Set<String> SV_prots = new HashSet<>();
    int min_skipped_exon = Integer.MAX_VALUE;
    int max_skipped_exon = Integer.MIN_VALUE;
    int min_skipped_bases = Integer.MAX_VALUE;
    int max_skipped_bases = Integer.MIN_VALUE;

    public ESSE(String gene_id, String gene_symbol, String chromosome, char strand, int svStart, int svEnd, int skipped_exon_count, int skipped_bases_count) {
        this.gene_id = gene_id;
        this.gene_symbol = gene_symbol;
        this.chromosome = chromosome;
        this.strand = strand;
        this.svStart = svStart;
        this.svEnd = svEnd;
        setSkippedExon(skipped_exon_count);
        setSkippedBases(skipped_bases_count);
    }
    public ESSE(Gene gene, Intron svIntron){
        this.gene = gene;
        this.gene_id = gene.id;
        this.gene_symbol = gene.symbol;
        this.chromosome = gene.chromosome;
        this.strand = svIntron.strand;
        this.svStart = svIntron.start;
        this.svEnd = svIntron.end;
    }
    public void setSkippedExon(int skipped_exon_count){
        if (skipped_exon_count > this.max_skipped_exon) this.max_skipped_exon = skipped_exon_count;
        if (skipped_exon_count < this.min_skipped_exon) this.min_skipped_exon = skipped_exon_count;
    }
    public void setSkippedBases(int skipped_bases_count){
        if (skipped_bases_count > this.max_skipped_bases) this.max_skipped_bases = skipped_bases_count;
        if (skipped_bases_count < this.min_skipped_bases) this.min_skipped_bases = skipped_bases_count;
    }
    public String getInterval(){
        return svStart + ":" + svEnd;
    }
    public String getWTIntrons(){
        return WT_introns.stream()
                .sorted(Comparator.comparing(c -> c.start))
                .map(c -> c.start + ":" + c.end)
                .reduce((a,b) -> a + "|" + b)
                .orElse("");
    }
    public void addAllWTIntron(List<Intron> wti){
        for (Intron i1 : wti) {
            boolean exists = false;
            for (Intron i2 : WT_introns) {
                if (i2.equals(i1)) {
                    exists = true;
                    break;
                }
            }
            if (!exists) WT_introns.add(i1);
        }
    }
    public void addAllWTExon(List<CDS> wtcds){
        for (CDS cds1 : wtcds) {
            boolean exists = false;
            for (CDS cds2 : WT_exons) {
                if (cds1.equals(cds2)) {
                    exists = true;
                    break;
                }
            }
            if (!exists) WT_exons.add(cds1);
        }
    }

    public String getProtsSV(){
        return String.join("|", SV_prots.stream().sorted().toList());
    }
    public String getProtsWT(){
        return String.join("|", WT_prots.stream().sorted().toList());
    }


}
