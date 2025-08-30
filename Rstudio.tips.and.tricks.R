# Rstudio tips and tricks

# ** Auto completion for functions **
# Below, typing "hea" should suggest the commonly used function "head()"
hea

# ** Auto completion for arguments and variable names **
# Typing out "palmerpeng" should suggest the package "palmerpenguins".
library(palmerpeng)

# ** Auto wrapping of quotes, parentheses, and `backticks` **
# Double clip the word blue below, then typing the quote key,
# just once will wtap the word in quotes.
plot(y = penguins$body_mass_g, x = penguins$species, col = blue)

# ** Comand-C to copy selection **
# ** Comand-X to cut (copy and delete) selection **
# ** Command-V to paste selection **
# ** Command-Z to undo an edit **

# ** Shift-Command-C to comment a line **
# Hit it mutliple times to toggle commenting on and off.
ls()

# ** Command-D to delete a line **

# ** Shift-Command-D to duplicate a line **
# Or duplicate a selection within or across lines!

# ** In Rstudio don't use Shift-Command-S for "Save-As"!
# Rstudio default is that Shift-Command-S is "Source"
# which runs all code in the active document -- usually not good! **

# ** Instead use Option-Command-S for "Save As" **

# Happy coding!
