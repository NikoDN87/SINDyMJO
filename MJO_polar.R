#########################################################################################
#------------------- Función para calcular Fase y Amplitud de MJO ----------------------#
#----------------------- a partir de componentes RMM1 y RMM2      ----------------------#
#                                   NDN - 25/02/22                                      #
#                                                                                       #
# Input: la función toma compo inputs las componentes x = RMM1 e y = RMM2               #
# Output: la función da como salida la fase discretizada (phi) y  amplitud (A)          #
#########################################################################################

MJO.polar = function(x,y) {
  
  A = sqrt(x^2 + y^2)        # amplitud de MJO
  df = data.frame(phi = atan2(y, x)*180/pi)         # ángulo entre (-180,180] de MJO
  df$phi = ifelse(df$phi<0, 360+df$phi, df$phi)     # paso a ángulo polar entre [0,360)
  # tengo que pasar phi a variable discreta con valore 1:8
  df$binned = cut(df$phi, breaks = seq(0,360,45), labels = c(seq(5,8),seq(1,4)))
  
  out.df = data.frame(Amplitude = A, phi = df$binned)  # armo data.frame de salida
  return(out.df)
  
}
