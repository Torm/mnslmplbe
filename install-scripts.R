# Build and install package.
# Make sure that target x86_64-pc-windows-gnu is installed. (rustup target install x86_64-pc-windows-gnu)
# If toolchain stable-msvc does not work, try stable-gnu. (rustup toolchain install stable-gnu)
# Try to restart RStudio after installing targets and toolchains.
devtools::document() # Build the package.
devtools::install() # Install the package locally.


# Other commands
devtools::load_all()
# Update rextendr and upgrade project
install.packages("rextendr")
rextendr::update_scaffold()
rextendr::use_crate("extendr-api", version = "0.9.0")
rextendr::vendor_crates()
rextendr::use_extendr()


library(rextendr)
rust_sitrep()
