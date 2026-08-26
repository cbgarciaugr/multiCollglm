# Artículos de reproducción (pkgdown)

Esta carpeta no contiene vignettes del paquete: son **artículos de
pkgdown**. Se renderizan solo para la web del paquete
(`pkgdown::build_site()` / `pkgdown::build_articles()`) y **no** se
compilan durante `R CMD build`/`R CMD check` ni se instalan con el
paquete — por eso `vignettes/articles` está en `.Rbuildignore`. Esto
permite que cada artículo use datasets externos o tarde más en ejecutarse
sin afectar al build/check del paquete ni a CRAN.

## Cómo añadir una nueva reproducción

1. Copia `reproduccion-plantilla.Rmd` con un nombre descriptivo, por
   ejemplo:

   ```
   reproduccion-perez-2019-consumo-energia.Rmd
   ```

2. Rellena las secciones siguiendo la plantilla (contexto, datos,
   diagnóstico original tal y como se hizo en el trabajo, diagnóstico con
   `multiCollglm`, comparación, conclusión, referencias).

3. Antes de incluir datos de terceros, comprueba su licencia. Si no se
   pueden redistribuir, usa una simulación fiel a sus propiedades
   estadísticas y dilo explícitamente en la sección "Datos".

4. Añade el artículo a `_pkgdown.yml` (sección `articles`, grupo
   "Reproducciones") para que aparezca listado en el menú de la web.

5. Genera la web localmente para revisar el resultado antes de publicar:

   ```r
   pkgdown::build_site()
   # o solo los artículos:
   pkgdown::build_articles()
   ```
