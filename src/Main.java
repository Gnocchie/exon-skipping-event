

import org.apache.commons.cli.*;

import java.io.IOException;


public class Main {
    public static void main(String[] args) throws IOException {



        Options options = new Options();
        // Updated to use long options with double dashes
        options.addOption(Option.builder()
                .longOpt("gtf")
                .hasArg()
                .desc("Path to gtf file")
                .build());
        options.addOption(Option.builder()
                .longOpt("o")
                .hasArg()
                .desc("Output path")
                .build());

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println("Error: " + e.getMessage());
            formatter.printHelp("java -jar exon_skipping.jar", options);
            return;
        }

        // Extract arguments using long options
        String gtf_path = cmd.getOptionValue("gtf");
        String out_path = cmd.getOptionValue("o");

        if (gtf_path == null || out_path == null) {
            System.out.println("Both --gtf and --o options are required.");
            formatter.printHelp("java -jar exonskipping.jar", options);
            return;
        }

        ESSE_Finder program = new ESSE_Finder();
        program.findEvents(gtf_path, out_path);

    }
}