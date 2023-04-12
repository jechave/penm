# from: https://emilkirkegaard.dk/en/2017/01/renaming-functions-in-r-packages-using-roxygen2/

#function generator
defunct = function(msg = "This function is depreciated") function(...) return(stop(msg))

#' @export
old_name = defunct("old_name changed name to new_name")



#' @export
delta_energy_dvm = defunct("delta_energy_dvm changed name to ddg_dv")

#' @export
delta_energy_dg_entropy = defunct("delta_energy_dg_entropy changed name to ddg_tds")
