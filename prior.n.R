prior.n <- function (n,equal=T) {

	if (!equal) {

		if (n>0.5 || (n>-0.5 && n<0)) {

			return(0.1-1/5/4)#assign lower probablity for speciation with highly asymetric split of niche

		} else {

			return(0.1+1/5/6)

		}

	} else {

		return(0.1)

	}

}
