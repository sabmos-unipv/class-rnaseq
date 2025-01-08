###################################
# Git clone delle necessarie repo #
###################################


git clone https://github.com/lescai-teaching/dataset_tutoring_rnaseq01.git # Dataset primo tutorato
git clone https://github.com/lescai-teaching/dataset_tutoring_rnaseq02.git # Dataset secondo tutorato


# Creazione della cartella di lavoro
cd /workspace/class-rnaseq
mkdir -p analysis_tutoring02
cd analysis_tutoring02


# Symbolic links alle reads originali
mkdir -p reads
cd reads
ln -s /workspace/class-rnaseq/dataset_tutoring_rnaseq02/raw_data/* .



###################
# Quantificazione #
###################

# Provare a eseguire salmon digitando "salmon" sul terminale
# A volte il comando potrebbe fallire sul terminale di RStudio su Codespaces
# Per risolvere digitare il comando: export PATH=${PATH}:/usr/local/bin


# Ciclo attraverso tutti i file che terminano con "_1.fasta.gz"
for sample in *_1.fasta.gz
do
    # Definizione dell'indice di riferimento
    index="/workspace/class-rnaseq/datasets_reference_only/trascriptome/chr21_transcripts_index"
    
    # Rimozione della parte finale "_1.fasta.gz" per ottenere il nome base
    name=${sample%_1.fasta.gz}
    
    # Messaggio di inizio elaborazione
    echo "Quantifying $name"
    
    # Salmon per quantificare l'espressione
    salmon quant \
        -p 2 \
        -i "$index" \
        -l IU \
        -1 "${name}_1.fasta.gz" \
        -2 "${name}_2.fasta.gz" \
        --validateMappings \
        -o "${name}.quant"
    
    # Messaggio di fine elaborazione
    echo -e "$name done now\n"
done



#########################################
# Esamina di un file di quantificazione #
#########################################


cd sample_01.quant
head quant.sf
