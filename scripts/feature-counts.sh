#!/bin/bash
# feature-counts.sh

script_name="feature-counts.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --File path to alignments, space seperated , control alignments followed by comparison groups"
  echo "-s  --Type of reads (e.g paired-end, single-end)"
  echo "-a  --Path to the annotation file"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.bam' -s 'single-end' -a '/path/to/gtf/gencode.v19.annotation.gtf' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}


main(){

    # Load required modules

    module load subread/1.6.3

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:a:o:vh opt
        do
            case $opt in
                f) alignments=$OPTARG;;
                s) sequencing_type=$OPTARG;;
                a) transcriptome=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $out ]] || [[ -z $alignments ]] || [[ -z $sequencing_type ]] || [[ -z $transcriptome ]]; then
        usage
    fi

    if [ $sequencing_type != 'single-end' ] && [ $sequencing_type != 'paired-end' ]; then
        usage
    fi

    if [ -z $out ]; then

        one_parent=$(dirname "${alignments}")
        out_dir=$(dirname "${one_parent}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi


    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

   # Define the output file name, based on the experiment
    #output_fp=${experiments}

    if [ ! -e $out_dir/featurecounts.txt ]; then

    if [ $sequencing_type = 'single-end' ]; then

    featureCounts -T 5 -t exon -g gene_id -a $transcriptome -o $out_dir/featurecounts.txt $alignments

    else

    featureCounts -p -T 5 -t exon -g gene_id -a $transcriptome -o $out_dir/featurecounts.txt $alignments

    fi
      # Get input and output files and then print out metadata.json file
        alignments_array=(${alignments//[,| ]/ })
        input_files=("$transcriptome" "${alignments_array[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_files=($out_dir\/featurecounts.txt $out_dir\/featurecounts.txt.summary)
        printf -v output "\"%s\"," "${output_files[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$output_fp\.metadata.json
        echo "* Finished."

    else
        echo "* $alignments has featurecounts generated"
    fi

    }

main "$@"
