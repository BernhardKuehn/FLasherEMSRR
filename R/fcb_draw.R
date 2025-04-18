# Classes and methods for plotting FCB diagrams
# fcb_draw.R
# FLasher/R/fcb_draw.R

# Copyright European Union, 2016
# Author: Finlay Scott (EC JRC) <finlay.scott@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

#' A base class for drawing FCB bits
#'
#' Look up and see the flying saucers cruising in the sky
#' I saw one myself it ain't no lie
#' Look down and see the road you're on as if you are on a marathon
#' That's the spirit, victory or die
#' 
#' @name basicBlock
#' @rdname basicBlock
#' @aliases basicBlock basicBlock-class
#'
#' @section Slots:
#'     \describe{
#'     \item{height}{height of box}
#'     \item{width}{width of box}
#'     \item{x_centre}{x coordinate of centre of box}
#'     \item{y_centre}{y coordinate of centre of box}
#'     \item{name}{text to go in the box}
#'     \item{name_cex}{size of text to go in the box}
#' }
#' @author Finlay Scott - EC JRC.
#' @keywords classes
setClass("basicBlock",
    representation(
        height = "numeric",
        width = "numeric",
        x_centre = "numeric",
        y_centre = "numeric",
        name = "character",
        name_cex = "numeric"
    ),
    prototype(
        height = NA_real_,
        width =  NA_real_,
        x_centre =  NA_real_,
        y_centre =  NA_real_,
        name = NA_character_,
        name_cex =  NA_real_
    )
)

#' A class for drawing a biological stock
#'
#' @name biolBlock
#' @rdname biolBlock
#' @aliases biolBlock biolBlock-class
#'
#' @section Slots:
#'     \describe{
#'     \item{height}{height of box}
#'     \item{width}{width of box}
#'     \item{x_centre}{x coordinate of centre of box}
#'     \item{y_centre}{y coordinate of centre of box}
#'     \item{name}{text to go in the box}
#'     \item{name_cex}{size of text to go in the box}
#'     \item{neck_length}{neck length}
#'     \item{circle_cex}{size of the circle bit}
#' }
#' @author Finlay Scott - EC JRC.
#' @keywords classes
setClass("biolBlock",
    representation(
        neck_length = "numeric",
        circle_cex = "numeric"
    ),
    prototype(
        neck_length =  NA_real_,
        circle_cex =  NA_real_
    ),
    contains = "basicBlock"
)

#' A class for drawing a fishery
#'
#' @name fisheryBlock
#' @rdname fisheryBlock
#' @aliases fisheryBlock fisheryBlock-class
#'
#' @section Slots:
#'     \describe{
#'     \item{height}{height of box}
#'     \item{width}{width of box}
#'     \item{x_centre}{x coordinate of centre of box}
#'     \item{y_centre}{y coordinate of centre of box}
#'     \item{name}{text to go in the box}
#'     \item{name_cex}{size of text to go in the box}
#'     \item{neck_length}{neck length}
#'     \item{no_tails}{number of catches}
#'     \item{tail_gap}{space between catches}
#'     \item{tail_length}{tail length}
#' }
#' @author Finlay Scott - EC JRC.
#' @keywords classes
setClass("fisheryBlock",
    representation(
        no_tails = "integer",
        neck_length = "numeric",
        tail_gap = "numeric",
        tail_length = "numeric"
    ),
    prototype(
        no_tails = NA_integer_,
        branch_length =  NA_real_,
        neck_length =  NA_real_,
        tail_gap = NA_real_,
        tail_length = NA_real_
    ),
    contains = "basicBlock"
)

#' A class for drawing a catch
#'
#' @name catchBlock
#' @rdname catchBlock
#' @aliases catchBlock catchBlock-class
#'
#' @section Slots:
#'     \describe{
#'     \item{height}{height of box}
#'     \item{width}{width of box}
#'     \item{x_centre}{x coordinate of centre of box}
#'     \item{y_centre}{y coordinate of centre of box}
#'     \item{name}{text to go in the box}
#'     \item{name_cex}{size of text to go in the box}
#'     \item{tail_length}{tail length}
#' }
#' @author Finlay Scott - EC JRC.
#' @keywords classes
setClass("catchBlock",
    representation(
        tail_length = "numeric"
    ),
    prototype(
        tail_length = NA_real_
    ),
    contains = "basicBlock"
)

#' A class for drawing the FCB matrix
#'
#' @name FCBDrawing
#' @rdname FCBDrawing
#' @aliases FCBDrawing FCBDrawing-class
#'
#' @section Slots:
#'     \describe{
#'     \item{biolBlocks}{A list of \link{biolBlock}s}
#'     \item{fisheryBlocks}{A list of \link{fisheryBlock}s}
#'     \item{catchBlocks}{A list of \link{catchBlock}s}
#'     \item{CBConnectors}{A list of coordinates connecting catches and biols}
#'     \item{FCB}{The FCB matrix}
#' }
#' @author Finlay Scott - EC JRC.
#' @keywords classes
setClass("FCBDrawing",
    representation(
        biolBlocks = "list",
        fisheryBlocks = "list",
        catchBlocks = "list",
        CBConnectors = "list",
        FCB = "matrix"),
    prototype(
        biolBlocks = list(),
        fisheryBlocks = list(),
        catchBlocks = list(),
        CBConnectors = list(),
        FCB = matrix()
    )
)

#' Generic method for drawing components the FCB matrix
#'
#' Used to draw individual components or the whole FCB matrix
#'
#' @param object The object - link{basicBlock}, link{biolBlock}, link{fisheryBlock}, link{catchBlock}, link{matrix}, link{FCBDrawing}
#' @param fisheryNames A vector of names of the FLFishery blocks
#' @param catchNames A vector of names of the FLCatch blocks
#' @param biolNames A vector of names of the FLBiol blocks
#' @param ... Other arguments
#' @return Nothing. Just draws something
#' @aliases draw draw-method
#' @rdname draw-FCB
setGeneric("draw", function(object, ...){
    standardGeneric("draw")
})

#' Method for drawing the basic block of the FCB matrix 
#'
#' @aliases draw-basicBlock
#' @rdname draw-FCB
setMethod("draw", signature(object="basicBlock"),
    function(object){
        # Draw the box round centre with height and width
        rect(object@x_centre - (object@width/2),
            object@y_centre - (object@height/2),
            object@x_centre + (object@width/2),
            object@y_centre + (object@height/2))
        # Add the text
        text(x=object@x_centre, y=object@y_centre, labels=object@name, cex=object@name_cex)
})

#' Method for drawing the biological block of the FCB matrix 
#'
#' @aliases draw-biolBlock
#' @rdname draw-FCB
setMethod("draw", signature(object="biolBlock"),
    function(object){
        # Call the basic drawing method
        callNextMethod()
        # Add the neck
        lines(x=c(object@x_centre, object@x_centre), y=c(object@y_centre + (object@height/2), object@y_centre + (object@height/2) + object@neck_length)) 
        # Add the circle 
        points(x=object@x_centre, y=object@y_centre + (object@height/2) + object@neck_length, pch=21, cex=object@circle_cex, col="black", bg="white")
})

#' Method for drawing the fishery block of the FCB matrix 
#'
#' @aliases draw-fisheryBlock
#' @rdname draw-FCB
setMethod("draw", signature(object="fisheryBlock"),
    function(object){
        # Call the basic drawing method
        callNextMethod()
        # Add the neck
        lines(x=c(object@x_centre, object@x_centre), y=c(object@y_centre - (object@height/2), object@y_centre - (object@height/2) - object@neck_length)) 
        # Add the branch
        branch_y <- object@y_centre - (object@height/2) - object@neck_length
        branch_length <- (object@no_tails -1) * object@tail_gap
        branch_x1 <- object@x_centre - branch_length/2
        branch_x2 <- object@x_centre + branch_length/2
        lines(x=c(branch_x1, branch_x2), y=c(branch_y, branch_y))
        # Or set tail gap as member instead of branch length?
        # Add the tails
        for (tail_count in 1:object@no_tails){
            tail_x <- branch_x1 + ((tail_count - 1) * object@tail_gap)
            lines(x=c(tail_x, tail_x), y=c(branch_y, branch_y - object@tail_length))
        }
})

#' Method for drawing the catch block of the FCB matrix 
#'
#' @aliases draw-catchBlock
#' @rdname draw-FCB
setMethod("draw", signature(object="catchBlock"),
    function(object){
        # Call the basic drawing method
        callNextMethod()
        # Add the tail
        lines(x=c(object@x_centre, object@x_centre), y=c(object@y_centre - (object@height/2), object@y_centre - (object@height/2) - object@tail_length)) 
})

#' Method for drawing the \link{FCBDrawing} object
#'
#' @aliases draw-FCBDrawing
#' @rdname draw-FCB
setMethod("draw", signature(object="FCBDrawing"),
    function(object){
        plot(x=c(0,1), y=c(0,1), type="n", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", ,bty="n")
        invisible(lapply(object@catchBlocks, draw))
        invisible(lapply(object@fisheryBlocks, draw))
        invisible(lapply(object@CBConnectors, function(coords){
            lines(x=c(coords["catch_x"], coords["biol_x"]), y = c(coords["catch_y"], coords["biol_y"]))
        }))
        invisible(lapply(object@biolBlocks, draw))
})

#' Method for drawing the FCB matric
#'
#' @aliases draw-matrix
#' @rdname draw-FCB
#' @examples
#' FCB <- matrix(c(1,1,1,1,2,2,2,1,2,2,2,3,2,2,4), nrow=5, byrow=TRUE)
#' draw(FCB)
setMethod("draw", signature(object="matrix"),
    function(object, fisheryNames=NULL, catchNames=NULL, biolNames=NULL){
        fcbd <- FCBDrawing(object, fisheryNames, catchNames, biolNames)
        draw(fcbd)
})

#' Method for drawing the FCB matric in the \link{fwdControl}
#'
#' @aliases draw-fwdControl
#' @rdname draw-FCB
setMethod("draw", signature(object="fwdControl"),
    function(object, fisheryNames=NULL, catchNames=NULL, biolNames=NULL){
        fcbd <- FCBDrawing(object@FCB, fisheryNames, catchNames, biolNames)
        draw(fcbd)
})


#' Generic constructor for making an \link{FCBDrawing} object
#'
#' @param FCB The (FCB \link{matrix})
#' @param fisheryNames A vector of names of the FLFishery blocks
#' @param catchNames A vector of names of the FLCatch blocks
#' @param biolNames A vector of names of the FLBiol blocks
#' @param ... Other arguments
#' @return An \link{FCBDrawing} object
#' @aliases FCBDrawing FCBDrawing-method
#' @rdname FCBDrawing
setGeneric("FCBDrawing", function(FCB, ...){
    standardGeneric("FCBDrawing")
})

#' Constructor for making an \link{FCBDrawing} object
#'
#' Has to figure out the height, width and positions of the FCB components
#' @aliases FCBDrawing-matrix
#' @rdname FCBDrawing
setMethod("FCBDrawing", signature(FCB="matrix"),
    function(FCB, fisheryNames=NULL, catchNames=NULL, biolNames=NULL){
        out <- new("FCBDrawing", FCB=FCB)
        # Sort out the FCB
        FCB <- data.frame(FCB)
        colnames(FCB) <- c("F","C","B")
        # FCB must be sorted by fishery and catch
        FCB <- FCB[order(FCB[,"F"],FCB[,"C"]),,drop=FALSE]
        # Need to add catch position in list - catches can catch multiple biols
        unique_FCB <- unique(FCB[,c("F","C")])
        unique_FCB <- cbind(unique_FCB,catch_posn = 1:nrow(unique_FCB))
        FCB <- merge(FCB, unique_FCB)
        # Note that FCB is now a data.frame

        # Ugh... loads of mucking about to get the right height, widths and lengths depending on FCB
        plot_width <- 1
        plot_height <- 1
        # Most number of blocks in a layer
        # Number of blocks in each layer
        no_fisheries <- max(FCB[,"F"])
        no_catches <- sum(tapply(FCB[,"C"],FCB[,"F"],max))
        no_biols <- max(FCB[,"B"])
        max_blocks <- max(no_fisheries, no_catches, no_biols)
        # Width based on biggest layer
        block_gap_width_ratio <- 2.0
        width <- 1 / (max_blocks + max_blocks / block_gap_width_ratio - 1)

        # Height cannot be too big, else cannot fit vertical space in
        height_neck <- 4 # The ideal height / neck ratio - ideally, height is 4 x neck
        # max height / neck ratio is 5/3
        # This then makes the height of each box (based on 3 * boxes + 5 * necks)
        max_height <- (plot_height * height_neck) / (3 * height_neck + 5) # height if neck = height / height_neck
        # If we use this box height, he neck length is
        max_neck_length <- (plot_height - 3 * max_height) / 5 # actual neck length to accomodate the height
        # This might be too small so we set a minimum neck length
        min_neck_length <- 0.1
        neck_length <- max(min_neck_length, max_neck_length)
        # Base height on the neck length
        height <- (plot_width - 5 * neck_length) / 3
        # This height might be too big for the width, conversely the width might be too big for the height
        # We take the min
        height <- min(height, width)
        width <- height
        tail_length <- neck_length

        # Fix names
        if(is.null(biolNames)){
            biolNames <- paste("Biology", 1:no_biols, sep="")
        }
        if(is.null(fisheryNames)){
            fisheryNames <- paste("Fishery", 1:no_fisheries, sep="")
        }
        if(is.null(catchNames)){
            catchNames <- paste("Catch", 1:no_catches, sep="")
        }

        # Make the catches
        catches <- list()
        gap_width <- (plot_width - (no_catches * width)) / (no_catches +1)
        y_centre <- height + (3 * neck_length) + (height / 2)
        # Spread boxes equally across plotting space
        x_centres <- ((1:no_catches) * gap_width) + (((1:no_catches)-1) * width) + (width / 2)
        for (i in 1:no_catches){
            name <- catchNames[i]
            nchar_name <-  nchar(name)
            name_cex <- (width * 40) / nchar_name
            catches[[i]] <- new("catchBlock", height = height, width = width, x_centre = x_centres[i], y_centre = y_centre, tail_length=tail_length, name = name, name_cex=name_cex) 
        }

        # Make the fisheries - based on catches position
        fisheries <- list()
        y_centre <- (2.5 * height) + (5 * neck_length)
        for (i in 1:no_fisheries){
            # x_centre depends on number of catches in fishery and total number of catches and position in list
            # Catches out of all catches that fishery has
            FCs <- unique(FCB[FCB[,"F"]==i,"catch_posn"])
            # Block should sit between them
            catches_xrange <- unlist(lapply(catches[FCs], function(x) x@x_centre))
            x_centre <- ((max(catches_xrange) - min(catches_xrange)) / 2) + catches_xrange[1]
            # tail_gap
            tail_gap <- (max(catches_xrange) - min(catches_xrange)) / (length(FCs)-1)
            if(is.nan(tail_gap)){ # If single catch
                tail_gap <- 0.0
            }
            name <- fisheryNames[i]
            nchar_name <-  nchar(name)
            name_cex <- (width * 40) / nchar_name
            fisheries[[i]] <- new("fisheryBlock", height = height, width = width, x_centre = x_centre, y_centre = y_centre, no_tails = length(FCs), tail_gap = tail_gap, neck_length=neck_length, tail_length=tail_length, name = name, name_cex=name_cex) 
        }

        # Make the Biols
        biols <- list()
        y_centre <- (height/2)
        circle_cex <- width * 20
        gap_width <- (plot_width  - (no_biols * width)) / (no_biols +1)
        # Spread boxes equally across plotting space
        x_centres <- ((1:no_biols) * gap_width) + (((1:no_biols)-1) * width) + (width / 2)
        # name_cex based on nchar
        for (i in 1:no_biols){
            name <- biolNames[i]
            nchar_name <-  nchar(name)
            name_cex <- (width * 40) / nchar_name
            biols[[i]] <- new("biolBlock", height = height, width = width, x_centre = x_centres[i], y_centre = y_centre, neck_length = neck_length, circle_cex = circle_cex, name = name, name_cex=name_cex) 
        }

        # Connections
        # List of coord arrays: x1,y1, x2,y2
        connectors <- list()
        for (i in 1:nrow(FCB)){
            # Position in list of biols - straightforward
            biol <- biols[[FCB[i,"B"]]]
            catch <- catches[[FCB[i,"catch_posn"]]]
            catch_y <- catch@y_centre - (catch@height/2) - catch@tail_length
            catch_x <- catch@x_centre
            biol_x <- biol@x_centre
            biol_y <- biol@y_centre + (biol@height/2) + biol@neck_length 
            connectors[[i]] <- c(catch_x = catch_x, catch_y = catch_y, biol_x = biol_x, biol_y = biol_y)
        }

        out@catchBlocks <- catches
        out@fisheryBlocks <- fisheries
        out@biolBlocks <- biols
        out@CBConnectors <- connectors
        return(out)
})


jukebox <- function(){
    tracks <- list(PinkFloyd_VegetableMan = c("In yellow shoes I get the blues", "Though I walk the streets with my plastic feet", "With my blue velvet trousers, make me feel pink", "There's a kind of stink about blue velvet trousers", "In my paisley shirt I look a jerk", "And my turquoise waistcoat is quite out of sight", "But oh oh my haircut looks so bad", "Vegetable man, where are you?", "So I've changed my dear, and I find my knees", "And I covered them up with the latest cut", "And my pants and socks all point in a box", "They don't make long of my nylon socks", "The watch, black watch", "My watch with a black face", "And a big pin, a little hole", "And all the lot is what I got", "It's what I wear, it's what you see", "It must be me, it's what I am", "Vegetable man", "I've been looking all over the place for a place for me", "But it ain't anywhere, it just ain't anywhere", "Vegetable man, vegetable man", "He's the kind of person, you just gonna see him if you can", "Vegetable man"),
                    OM_Gethsemane = c("Sing the advatic - penitent rite to send","Toward redemption - strength attends the seeker", "A sovereign dreadnaught breathes and approach the state serene", "Interiorizer of the mind from world labyrinth", "Nocodemus awaits in vigil weeping", "The Arahat rising and the healing ghost descends", "Lamentations cease enter rarefied light prevails", "Devekut gleams sing freedom from tamasic field", "And Ezekiel saw the wheel", "Prana exits toward the azure locus", "On vrittiless path he walks", "Contemplative ascends", "Temple-cave of the heart shrine", "Vigil on the mountain"),
                    Clutch_BeholdTheColossus = c("Drive us by pitchfork to caves by the sea", "Atop the mountain, under cities", "We are monsters", "We are older than men", "Yes, we are Gorgons", "And we will have our revenge", "Frighten their children by speaking our name", "Epoch and eon, it's always the same", "We are monsters", "But we will no longer run", "We are the Cyclops", "And our time has come", "This world was ours before it had ever known men", "Soon it will be ours again", "Grant them their tales of their forefather's glory", "For soon it will be ours again", "Atop the buttress the gargoyle sits", "The fragile foundation crumbles to bits", "We are monsters and they are but men", "Behold the Colossus rises again!"),
                    OM_StateOfNonReturn = c("Traveler now reach the stream. The astral flight adapter", "From the pain-sheath life ascends - the Non-returner sees", "Empathy release me - and the phoenix rise triumphant", "And walks onto the certitude ground - the soul's submergence ends", "From the rounds of rebirth - he arrives onto the deathless", "Light bores through the adjunct worlds - the soul-galleon prevails", "Liberates in wisdom to complete state of negation", "The five roads subsumed by grace - emancipates from dream"),
                    SatansSatyrs_CurseOfTheCorpse = c("When the unmanned belfry resounds with grim knell", "A fiend returns from halfway to Hell", "With human shape, but demon's gait", "And five grave-bound centuries of hunger to sate", "Do you not remember", "Or do you wish to forget?", "Tonight the corpse will rise", "This night of Walpurgis", "In life he practiced blood-sacrifice", "Left not untapped the blackest vice", "His occult prowess preserved his flesh", "Not among the living nor of the dead", "They say the dead hate the living", "But is it really true?", "You won't guess twice", "When they're coming after you", "You'll find no life behind those eyes", "Don't even risk a glance", "A stake through the heart is not enough", "Don't even stand a chance"),
                    OrangeGoblin_RedTideRising = c("Chaos is calling, the kingdoms are falling", "The planets collide in the sky", "The old ones command, by a curse or by hand", "We've come to the end of our times", "Thoughts paralyzed, a red tide will rise", "Dystopia seen through cephalopod eyes", "Creatures reborn, from weird tales untold", "Visions of horror, dimensions unknown to this world!", "Gods are descending, a force unrelenting", "Cthulhu will rise from the sea", "Turning blood red and awaking the dead", "From their aeons of suffering sleep", "Monoliths rise, to fiery skies", "Empires dead in the shadow of flies", "Serpents unwind, to torture and bind", "No one survives the ultimate fall of mankind!", "Atlantis will rise!", "The Kraken alive!", "These are Fortean times!", "No sanctuary", "No Hail Mary", "No visionary", "Can save you from what has begun!", "Infidels fall to their knees and they start", "Praying to something that's never been there in their hearts", "Searing white light from the guardians flesh", "Azrael rising, wings of the angel of death", "No sanctuary", "No Hail Mary", "No visionary", "Can save you from what has begun!", "This is Leviathan!", "Death of the sun!", "Age of darkness has won"),
                    Cud_PurpleLoveBalloon = c("The first time that I saw you", "You were running scared", "You assured me", "Your clothes were making chase", "Your mouth was making shapes before me", "The second time you touched me", "Your hands sent ripples through me", "You're one in a million, take me home", "Make dragons fly, make dragons roam", "You're one in a million", "And I will see you soon", "I want to take you high in my...", "Purple Love Balloon", "The second time you touched me", "Your hands went running through me", "It was white so I could find me", "The second time you loved me", "The dragon was inside me", "I want to take you high", "In my Purple Love Balloon", "I want to thumb a ride", "In my Purple Love Balloon", "I want to make you fly", "In my Purple Love Balloon"),
                    SuperFurryAnimals_DacwHi = c("Dacw hi, a'I gwyneb mewn poen", "Mae hi newydd gael ei phigo gan byr ar ei chroen", "Mae'n gafael mewn papur cyfagos a'I daro ar y wal", "Mae'n gadael gwennynen ar ?ydd yn sicr wedi ei dal", "Mae hi'n ddoniol", "Mae'n sylwi pob dim", "Ac mae'n medru gweld", "Drwy gefn ei phen", "Dim dioloch I'r drefn", "Mae'n fore iau ac mae'n mynd I'r gwaith", "Wrth gau ei drws mae hi'n tynnu at saith", "Wrth gau ei ch?hag yr oerni mae hi'n colli pishyn punt", "Mae'n ei bigo o fyny yn syth cyn troi n? mewn I'r gwynt", "Mae hi'n ddoniol", "Mae'n sylwi pob dim", "Ac mae'n medru gweld", "Drwy gefn ei phen", "Dim dioloch I'r drefn", "Mae hi'n ddoniol", "Mae'n sylwi pob dim", "Mae hi'n ddoniol", "Neith hi chwerthin ar ddim", "Mae hi'n ddoniol", "Mae'n sylwi pob dim", "Ac mae'n medru gweld", "Drwy gefn ei phen", "Dim dioloch I'r drefn", "Unwaith drachefn", "Amen"),
                    MonsterMagnet_ILiveBehindTheClouds = c("I've got a feeling that no one cares", "'Bout all that confetti I throw in the air", "Nothing's important yet everything is", "If there ain't no photo I just don't exist", "And the winds blow", "And the sky looks cool", "So I make my home in the clouds", "I see the oceans, the mountains, the plains", "The rats and the suckers, the strong and the lame", "Sometimes I dream about flying on down", "Falling in love and life on the ground", "And the winds blow", "And the sun comes up", "And I stay behind the clouds", "Yes I live behind the clouds", "Raise Surya to burn off the soil", "Make it safe for my eyes", "I'd pay good money just to get proved wrong", "For once in my goddamned life", "Pull me in, drag me down", "Into your ocean of love", "My wings are ragged and my soul is gone", "So tired of floating above", "Sing to me baby your sweet refrain", "I feel I'm coming apart in this rain", "That's all I need if that's all I hear", "Some beautiful music as I disappear", "And the winds blow", "And the sun goes down", "And I stay behind the clouds", "Yes I stay behind the clouds", "Yeah I live behind the clouds"))
    trackno <- round(runif(1,min=1,max=length(tracks)))
    lyricno <- round(runif(1,min=1,max=length(tracks[[trackno]])))
    return(tracks[[trackno]][lyricno])
}

.onAttach <- function(libname, pkgname ) {
    startup <- paste("FLasherEMSRR: ", jukebox(), sep="")
    packageStartupMessage(startup)
}
