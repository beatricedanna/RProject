# ==================================================================
# Dockerfile Semplice e Veloce per Progetto R Multi-omico
# Base: Bioconductor (ha già R, dipendenze di sistema e BiocManager)
# ==================================================================

# 1. Partiamo da un'immagine ufficiale di Bioconductor.
#    Questa immagine include già una versione recente di R e molte librerie di sistema.
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# 2. Imposta la directory di lavoro
WORKDIR /project

# 3. Installa le dipendenze R (inclusi Seurat e Signac)
#    Copiamo solo il file DESCRIPTION per ottimizzare la cache di Docker.
COPY DESCRIPTION .

# Definiamo un argomento per passare il token di GitHub in modo sicuro.
ARG GITHUB_PAT

# Installiamo TUTTE le dipendenze (da CRAN, Bioc, GitHub) in un colpo solo
# leggendo il file DESCRIPTION. Questo è il cuore della build.
# La base Bioconductor ha già 'remotes' pre-installato.
RUN --mount=type=secret,id=GITHUB_PAT \
  R -e 'Sys.setenv(GITHUB_PAT = scan("/run/secrets/GITHUB_PAT", what = "character", quiet = TRUE)); \
        BiocManager::install(ask = FALSE); \
        remotes::install_deps(dependencies = TRUE)'

# 4. Installa Git LFS (uno dei pochi comandi di sistema ancora necessari)
RUN apt-get update && apt-get install -y --no-install-recommends git-lfs && rm -rf /var/lib/apt/lists/*

# 5. Copia tutto il resto del progetto (codice, dati, etc.)
COPY . .

# 6. Scarica i dati pesanti da Git LFS
RUN git lfs install --force && git lfs pull

# 7. Installa il tuo pacchetto locale
# Il punto '.' significa "installa il pacchetto in questa directory".
RUN R CMD INSTALL .

# 8. (Opzionale ma consigliato) Pre-calcola la vignetta
RUN R -e "devtools::build_vignettes()"

# 9. Imposta la directory finale e il comando di default
WORKDIR /analysis
CMD ["R"]
