import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.SimpleTimeZone;

public class ReaderGTF {
    public static void read(String file, HashMap<String, Gene> genes) throws IOException {
        FileReader fr = new FileReader(file);
        BufferedReader br = new BufferedReader(fr);
        String line;
        while ((line = br.readLine()) != null) {

            if (line.startsWith("#")) continue;
            String[] columns = line.split("\t");

            String chromosome = columns[0];
            int start = Integer.parseInt(columns[3]);
            int end = Integer.parseInt(columns[4]);
            char strand = columns[6].charAt(0);
            String[] attributes = columns[8].split(";");

            String geneId = fetchAttribute("gene_id", attributes);
            String geneSymbol = fetchAttribute("gene_name", attributes);
            String transcriptId = fetchAttribute("transcript_id", attributes);
            if (geneId == null || transcriptId == null) continue;

            Gene gene = genes.computeIfAbsent(geneId, k -> new Gene(geneId,geneSymbol, strand, chromosome));
            Transcript transcript = gene.transcriptHashMap.computeIfAbsent(transcriptId, k -> new Transcript(transcriptId, gene));

            if (!columns[2].equals("CDS")) continue;

            String proteinId = fetchAttribute("protein_id", attributes);
            gene.proteinHashSet.add(proteinId);
            transcript.proteinList.add(proteinId);

            CDS cds = new CDS(transcript, proteinId, start, end, strand);
            transcript.cdsList.add(cds);
            gene.addUniqueCDS(cds, transcript);
        }
    }
    public static String fetchAttribute(String key, String[] attributes) {
        for (String attribute : attributes) {
            attribute = attribute.trim();
            if (attribute.startsWith(key)) {
                return attribute.split(" ")[1].replace("\"", "");
            }
        }
        return null;
    }
}
