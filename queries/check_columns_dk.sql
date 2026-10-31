SELECT specieskey, taxonomicstatus, GBIF_ISEA3HCode(6, decimallatitude, decimallongitude, 0) AS cell6, GBIF_ISEA3HCode(7, decimallatitude, decimallongitude, 0) AS cell7, CASE WHEN level0gid IS NULL THEN 'sea' ELSE 'land' END AS part, COUNT(*) AS n
FROM occurrence
WHERE countrycode = 'DK' AND hascoordinate = TRUE AND hasgeospatialissues = FALSE AND occurrencestatus = 'PRESENT'
GROUP BY specieskey, taxonomicstatus, GBIF_ISEA3HCode(6, decimallatitude, decimallongitude, 0), GBIF_ISEA3HCode(7, decimallatitude, decimallongitude, 0), CASE WHEN level0gid IS NULL THEN 'sea' ELSE 'land' END
