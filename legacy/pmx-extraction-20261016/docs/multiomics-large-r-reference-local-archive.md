# MultiOmics — résultats R volumineux et confidentialité

**Pour les chercheurs / For researchers.** La plateforme peut retourner les tableaux DESeq2 / limma-voom complets dans le dossier ZIP **tant que leur taille reste dans la limite JSON du backend R** (64 Mio par méthode dans la PR #81). Quand une sortie dépasse cette limite, MultiOmics ne doit **jamais la présenter comme exhaustive dans le ZIP**.

## Option 1 — fichiers ordinaires

Aucune configuration supplémentaire : l'API retourne `referenceTableComplete=true`, `referenceTableStatus=complete`, et le ZIP contient `reference_r/<méthode>.csv` avec **toutes** les lignes ajustées. Les `tableaux/` sont distincts : ce sont les résultats du moteur exploratoire navigateur, pas des résultats R.

## Option 2 — fichiers trop volumineux (stockage LOCAL facultatif)

**Avant de lancer le backend R**, choisissez sur l'ordinateur du chercheur un répertoire dédié, privé et doté d'un espace disque suffisant. Le site GitHub Pages ne peut pas accéder automatiquement à ce répertoire et ne doit pas téléverser son contenu ailleurs.

Sous Windows PowerShell :

```powershell
$archive = Join-Path $env:USERPROFILE "Documents\PMxMultiOmics\resultats-R"
New-Item -ItemType Directory -Force -Path $archive | Out-Null
$env:PMX_MULTIOMICS_RESULTS_DIR = $archive
Rscript multiomics-engine/run_backend.R
```

Sous Linux / macOS :

```bash
install -d -m 700 "$HOME/PMxMultiOmics/resultats-R"
PMX_MULTIOMICS_RESULTS_DIR="$HOME/PMxMultiOmics/resultats-R" \
  Rscript multiomics-engine/run_backend.R
```

Si le CSV R excède la taille JSON :

- L'API annonce **`referenceTableComplete=false`** et **`referenceTableStatus=too_large_archived_locally`** : il ne prétend pas être dans le ZIP.
- Un dossier aléatoire `pmx_reference_.../fitted_reference.csv` est créé dans le répertoire local configuré. Le backend conserve le nom du dossier (`referenceTableLocalArchiveId`) et l'empreinte SHA-256 (`referenceTableLocalArchiveSha256`), **pas le chemin absolu** dans sa réponse.
- Le dossier est créé avec les permissions 0700 et le CSV 0600 sur les systèmes POSIX. **Sous Windows, vérifier les ACL NTFS** du répertoire personnel : les modes POSIX ne remplacent pas des permissions Windows correctes.
- Le chercheur retrouve le CSV localement, puis vérifie son empreinte :

```powershell
Get-FileHash -Algorithm SHA256 "$archive\pmx_reference_...\fitted_reference.csv"
```

ou :

```bash
sha256sum "$HOME/PMxMultiOmics/resultats-R/pmx_reference_.../fitted_reference.csv"
```

- Le plafond supplémentaire de conservation par méthode est actuellement **1 Gio** ; au-delà, le moteur retourne une sortie incomplète **non conservée**. Aucun stockage par défaut, ni export vers une plateforme distante, n'est activé.
- **Si l'archivage n'est pas configuré ou échoue :** `referenceTableStatus=too_large`, CSV non retourné, *aucun* fichier à récupérer de l'API. La marche à suivre est de configurer un répertoire local et relancer. L'ancien message affirmant qu'un fichier subsistait après nettoyage des dossiers temporaires était trompeur : il a été corrigé.

## Recommandations de reproductibilité

Conserver : (1) le projet et les données originales, (2) les résultats R exhaustifs, (3) le code et la version des packages, (4) le contraste exact et les variables de confusion, (5) l'empreinte SHA-256 des fichiers. Une sortie du moteur R local peut contenir des informations médicales ou des identifiants internes ; vérifier le droit de partage avant toute transmission.

La présence d'une table complète et son SHA-256 sont des **preuves d'exhaustivité et d'intégrité technique**. Elles ne prouvent ni la qualité du protocole, ni le contrôle du FDR, ni la signification biologique.

*English summary:* Local retention is opt-in using `PMX_MULTIOMICS_RESULTS_DIR`. A file too large for the JSON response remains marked **incomplete in the browser/ZIP**; when locally archived, its random folder ID and SHA-256 allow retrieval and verification without transferring confidential data. No archive is created by default. Inspect Windows ACLs when applicable.
