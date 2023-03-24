# Import of the bibliography

  insert_head()
  
# container -------
  
  bib <- list()
  
# reading the BibTex ------
  
  insert_msg('Reading the BibTex')
  
  bib <- read_bib('./paper/markdown/pah_biblio.bib')
  
# END ------
  
  insert_tail()