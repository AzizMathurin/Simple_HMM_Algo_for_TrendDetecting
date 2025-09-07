#########Simple HMM for Detecting Bullish and Bearish Markets#####


##########forward_algorithm


forward_algorithm <- function(A, B, pi, O) {
  
  # Arguments :
  #   A  : Matrice de transition des états cachés (2x2)
  #   B  : Matrice de probabilité d'émission des observations (2x2)
  #   pi : Vecteur de probabilités initiales (2x1)
  #   O  : Vecteur des observations (caractères "u" ou "d")
  #
  # Retour :
  #   Une liste contenant :
  #   - probability : Probabilité de la séquence d'observations
  #   - alpha       : Matrice des probabilités partielles α
  
  # Exemple d'utilisation :
  # Définition des matrices de transition, d'émission et des probabilités initiales
  # A <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
  # B <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
  # pi <- c(0.6, 0.4)
  # O <- c("u", "d", "u", "u", "d")
  # Exécution de l'algorithme
  # result <- forward_algorithm(A, B, pi, O)
  
  # Vérification des entrées
  if (!is.matrix(A) || !is.matrix(B) || !is.vector(pi) || !is.vector(O)) {
    stop("Les entrées doivent être des matrices ou des vecteurs.")
  }
  if (nrow(A) != 2 || ncol(A) != 2) {
    stop("La matrice A doit être de dimension 2x2.")
  }
  if (nrow(B) != 2 || ncol(B) != 2) {
    stop("La matrice B doit être de dimension 2x2.")
  }
  if (length(pi) != 2) {
    stop("Le vecteur pi doit être de longueur 2.")
  }
  if (!all(O %in% c("u", "d"))) {
    stop("Les observations doivent être 'u' ou 'd'.")
  }
  # Nombre d'observations
  T <- length(O)
  
  # Initialisation de la matrice alpha (2 x T)
  alpha <- matrix(0, nrow = 2, ncol = T)
  
  # Initialisation de la première colonne de alpha
  if (O[1] == "u") {
    alpha[, 1] <- pi * B[, 1]  # Multiplication élément par élément
  } else  {
    alpha[, 1] <- pi * B[, 2]
  } 
  
  # Récursion pour remplir la matrice alpha
  for (t in 2:T) {
    alpha[, t] <- c(
      alpha[1, t - 1] * A[1, 1] + alpha[2, t - 1] * A[2, 1],
      alpha[1, t - 1] * A[1, 2] + alpha[2, t - 1] * A[2, 2]
    )
    
    if (O[t] == "u") {
      alpha[, t] <- alpha[, t] * B[, 1]  # Multiplication élément par élément
    } else {
      alpha[, t] <- alpha[, t] * B[, 2]
    } 
  }
  
  # Calcul de la probabilité de la séquence d'observations
  proba_O <- sum(alpha[, T])
  
  # Retourne la probabilité et la matrice alpha
  return(list(probability = proba_O, alpha = alpha))
}


#########Backward_algorithm


backward_algorithm <- function(A, B, pi, O) {
  
  # Arguments :
  #   A  : Matrice de transition des états cachés (2x2)
  #   B  : Matrice de probabilité d'émission des observations (2x2)
  #   pi : Vecteur de probabilités initiales (2x1)
  #   O  : Vecteur des observations (caractères "u" ou "d")
  #
  # Retour :
  #   Une liste contenant :
  #   - probability : Probabilité de la séquence d'observations
  #   - beta       : La matrice beta calculée par l'algorithme backward.
  
  # Exemple d'utilisation :
  # Définition des matrices de transition, d'émission et des probabilités initiales
  # A <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
  # B <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
  # pi <- c(0.6, 0.4)
  # O <- c("u", "d", "u", "u", "d")
  # Exécution de l'algorithme
  # result <- backward_algorithm(A, B, pi, O)
  
  # Vérification des entrées
  if (!is.matrix(A) || !is.matrix(B) || !is.vector(pi) || !is.vector(O)) {
    stop("Les entrées doivent être des matrices ou des vecteurs.")
  }
  if (nrow(A) != 2 || ncol(A) != 2) {
    stop("La matrice A doit être de dimension 2x2.")
  }
  if (nrow(B) != 2 || ncol(B) != 2) {
    stop("La matrice B doit être de dimension 2x2.")
  }
  if (length(pi) != 2) {
    stop("Le vecteur pi doit être de longueur 2.")
  }
  if (!all(O %in% c("u", "d"))) {
    stop("Les observations doivent être 'u' ou 'd'.")
  }
  
  # Nombre d'observations
  num_observations <- length(O)
  
  # Création de la matrice beta (2 x num_observations)
  beta <- matrix(0, nrow = 2, ncol = num_observations)
  
  # Initialisation de la dernière colonne de beta
  beta[, num_observations] <- 1
  
  # Récursion pour remplir la matrice beta
  for (t in seq(num_observations - 1, 1, by = -1)) {
    if (O[t + 1] == "u") {
      beta[1, t] <- sum(beta[, t + 1] * A[1, ] * B[, 1])
      beta[2, t] <- sum(beta[, t + 1] * A[2, ] * B[, 1])
    } else  {
      beta[1, t] <- sum(beta[, t + 1] * A[1, ] * B[, 2])
      beta[2, t] <- sum(beta[, t + 1] * A[2, ] * B[, 2])
    } 
  }
  
  # Calcul de la probabilité de la séquence d'observations
  if (O[1] == "u") {
    proba_O <- sum(beta[, 1] * pi * B[, 1])
  } else {
    proba_O <- sum(beta[, 1] * pi * B[, 2])
  }
  
  # Retourne la probabilité et la matrice beta
  return(list(probability = proba_O, beta = beta))
}


#####Decoding : Viterbi




viterbi_algorithm <- function(A, B, pi, O) {
  
  # Arguments :
  #   A  : Matrice de transition des états cachés (2x2)
  #   B  : Matrice de probabilité d'émission des observations (2x2)
  #   pi : Vecteur de probabilités initiales (2x1)
  #   O  : Vecteur des observations (caractères "u" ou "d")
  #
  # Retour :
  #   Une liste contenant :
  #   - chemin : Le chemin optimal des états cachés
  #   - probaOpti      : La probabilité du chemin optimal.
  
  # Exemple d'utilisation :
  # Définition des matrices de transition, d'émission et des probabilités initiales
  # A <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
  # B <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
  # pi <- c(0.6, 0.4)
  # O <- c("u", "d", "u", "u", "d")
  # Exécution de l'algorithme
  # result <- viterbi_algorithm(A, B, pi, O)
  
  # Vérification des entrées
  if (!is.matrix(A) || !is.matrix(B) || !is.vector(pi) || !is.vector(O)) {
    stop("Les entrées doivent être des matrices ou des vecteurs.")
  }
  if (nrow(A) != 2 || ncol(A) != 2) {
    stop("La matrice A doit être de dimension 2x2.")
  }
  if (nrow(B) != 2 || ncol(B) != 2) {
    stop("La matrice B doit être de dimension 2x2.")
  }
  if (length(pi) != 2) {
    stop("Le vecteur pi doit être de longueur 2.")
  }
  if (!all(O %in% c("u", "d"))) {
    stop("Les observations doivent être 'u' ou 'd'.")
  }
  
  # Longueur de la séquence d'observation
  num_observations <- length(O)
  num_states <- 2  # Nombre d'états
  
  # Création des matrices viterbi et Psi (2 x num_observations)
  viterbi <- matrix(0, nrow = num_states, ncol = num_observations)
  Psi <- matrix(0, nrow = num_states, ncol = num_observations)
  
  # Initialisation
  if (O[1] == "u") {
    viterbi[, 1] <- pi * B[, 1]  # Multiplication élément par élément
  } else {
    viterbi[, 1] <- pi * B[, 2]
  } 
  
  # Récursion
  for (t in 2:num_observations) {
    for (s in 1:num_states) {
      trans_probs <- viterbi[, t - 1] * A[, s]  # Probabilités de transition
      max_prob <- max(trans_probs)  # Meilleure probabilité
      max_state <- which.max(trans_probs)  # État optimal précédent
      
      viterbi[s, t] <- max_prob
      Psi[s, t] <- max_state
    }
    
    # Multiplication par la probabilité d'observation
    if (O[t] == "u") {
      viterbi[, t] <- viterbi[, t] * B[, 1]
    } else  {
      viterbi[, t] <- viterbi[, t] * B[, 2]
    } 
  }
  
  # Terminaison
  probaOpti <- max(viterbi[, num_observations])
  last_state <- which.max(viterbi[, num_observations])
  
  # Rétropropagation pour trouver le chemin optimal
  cheminOpti <- numeric(num_observations)
  cheminOpti[num_observations] <- last_state
  
  for (t in (num_observations - 1):1) {
    cheminOpti[t] <- Psi[cheminOpti[t + 1], t + 1]
  }
  
  # Retourne le chemin optimal et la probabilité du chemin optimal
  return(list(chemin = cheminOpti, probaOpti = probaOpti))
}


####Baum-Welch Forward-Backward algorithm


baum_welch_algorithm <- function(Alpha, Beta, O, P_O, A, B) {
  
  
  # Arguments :
  #   A  : Matrice de transition des états cachés (2x2)
  #   B  : Matrice de probabilité d'émission des observations (2x2)
  #   pi : Vecteur de probabilités initiales (2x1)
  #   O  : Vecteur des observations (caractères "u" ou "d")
  #   Alpha Matrice alpha calculée par l'algorithme forward.
  #   Beta Matrice beta calculée par l'algorithme backward.
  #   P_O Probabilité de la séquence d'observations.
  
  #  Alpha <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, byrow = TRUE)
  #  Beta <- matrix(c(0.4, 0.3, 0.2, 0.1), nrow = 2, byrow = TRUE)
  #  O <- c("u", "d", "u")
  #  P_O <- 0.05
  #  A <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
  #  B <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
  #  baum_welch_algorithm(Alpha, Beta, O, P_O_l, A, B)
  
  # Vérification des entrées
  if (!is.matrix(A) || !is.matrix(B) || !is.vector(pi) || !is.vector(O)) {
    stop("Les entrées doivent être des matrices ou des vecteurs.")
  }
  if (nrow(A) != 2 || ncol(A) != 2) {
    stop("La matrice A doit être de dimension 2x2.")
  }
  if (nrow(B) != 2 || ncol(B) != 2) {
    stop("La matrice B doit être de dimension 2x2.")
  }
  if (length(pi) != 2) {
    stop("Le vecteur pi doit être de longueur 2.")
  }
  if (!all(O %in% c("u", "d"))) {
    stop("Les observations doivent être 'u' ou 'd'.")
  }
  
  # Calcul de la matrice gamma
  gamma <- (Alpha * Beta) / P_O
  
  # Nombre d'états et longueur de la séquence d'observations
  number_states <- 2
  n <- length(O)
  
  # Initialisation de la matrice Epsilon
  Epsilon <- matrix(0, nrow = number_states * number_states, ncol = n - 1)
  
  # Calcul de la matrice Epsilon
  for (t in 1:(n - 1)) {
    compteur <- 0
    for (i in 1:number_states) {
      for (j in 1:number_states) {
        compteur <- compteur + 1
        Epsilon[compteur, t] <- Alpha[i, t] * A[i, j] *
          ifelse(O[t + 1] == "u", B[j, 1], B[j, 2]) * Beta[j, t + 1]
      }
    }
  }
  Epsilon <- Epsilon / P_O
  
  # Mise à jour de la matrice A
  vecteur_aij <- c(0, 0, 0, 0)
  a1 <- Epsilon[1:2, ]
  vecteur_aij[1:2] <- rowSums(a1) / sum(a1)
  a2 <- Epsilon[3:4, ]
  vecteur_aij[3:4] <- rowSums(a2) / sum(a2)
  A <- matrix(vecteur_aij, nrow = number_states, ncol = number_states, byrow = TRUE)
  
  # Mise à jour de la matrice B
  vectorH <- (O == "u")
  b1_H <- sum(gamma[1, vectorH]) / sum(gamma[1, ])
  b1_T <- sum(gamma[1, !vectorH]) / sum(gamma[1, ])
  b2_H <- sum(gamma[2, vectorH]) / sum(gamma[2, ])
  b2_T <- sum(gamma[2, !vectorH]) / sum(gamma[2, ])
  B <- matrix(c(b1_H, b2_H, b1_T, b2_T), nrow = 2, byrow = FALSE)
  
  # Mise à jour du vecteur pi
  pi <- gamma[, 1]
  
  # Retourne les matrices A, B et le vecteur pi mis à jour
  return(list(A = A, B = B, pi = pi))
}




######Main function


HHMM_2State <- function(A, B, pi, O,n=200){
  #n : le nombre d'iteration
  for (i in 1:n) {
    result1 <- forward_algorithm(A, B, pi, O)
    alpha <- result1$alpha
    Proba_O <- result1$probability
    
    result2 <- backward_algorithm(A, B, pi, O)
    beta <- result2$beta
    
    result3 <- baum_welch_algorithm(alpha, beta, O, Proba_O, A, B)
    A <- result3$A
    B <- result3$B
    pi <- result3$pi
  }
  
  result4 <- viterbi_algorithm(A, B, pi, O)
  chemin <- result4$chemin 
  
  return(list(A = A, B = B, pi = pi, chemin = chemin))
}

#####Application 1 Coin tossing game

Coin_flip=c('u','u','d','u','d',
            'd','u','d','u','u',
            'd','u','u','u','u',
            'd','d','u','d','d',
            'd','d','d','u','d',
            'd','u','d','u','d',
            'u','d','u','d','u',
            'u','d','u','u','u',
            'u','u','d','d','d',
            'd','u','d','u','d',
            'd','d','d','d','u',
            'd','u','d','d','u',
            'u','d','u','d','u',
            'u','u','u','d','u',
            'u','u','u','u','d',
            'u','d','u','u','d',
            'd','u','d','d','d',
            'u','d','d','d','u',
            'u','u','d','u','u',
            'd','u','u','u','d',
            'u','d','d','u','u',
            'u','u','u','u','u',
            'u','d','d','u','d',
            'd','d','u','d','u',
            'd','d','d','u','u',
            'u','d','d','u','d',
            'd','u','u','u','d')

realState=c(rep(1,15),rep(2,15),rep(1,10),rep(2,20),rep(1,20),rep(2,10),rep(1,20),rep(2,20),rep(1,5))

#Initialize matrices
A <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE)#initialisation A
B <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE)#initialisation B
pi <- c(0.5, 0.5)#initialisation pi




gamme1=HHMM_2State(A, B, pi, Coin_flip, 400)
abs(gamme1$chemin-realState)
accuracy=(1-sum(abs(gamme1$chemin-realState))/135)*100






# Initialisation des matrices by some guess
A <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)#initialisation A
B <- matrix(c(0.7, 0.3, 0.3, 0.7), nrow = 2, byrow = TRUE)#initialisation B
pi <- c(0.5, 0.5)#initialisation pi


gamme2=HHMM_2State(A, B, pi, Coin_flip, 400)
abs(gamme2$chemin-realState)
accuracy=(1-sum(abs(gamme2$chemin-realState))/135)*100
accuracy 


#####Detecting bearish and bullish markets


library("quantmod")


# Télécharger les données du Dow Jones Industrial Average (DJI)
DJI <- getSymbols(
  Symbols = "^DJI",
  src = "yahoo",
  from = "2023-01-01",
  to = "2024-01-01",
  periodicity = "daily",
  auto.assign = FALSE
)


# Extraire les prix ajustés
DJIPRICES <- DJI[, "DJI.Adjusted"]


# Calcul des moyennes mobiles
MA10 <- SMA(DJIPRICES, n = 10)
MA20 <- SMA(DJIPRICES, n = 20)


# Tracer le graphique des prix avec les moyennes mobiles
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))  # Diviser la fenêtre graphique en 2 parties


plot(DJIPRICES, main = "DJI Prices with Moving Averages", ylab = "Prices", xlab = "Date", col = "black")
lines(MA10, col = "blue", lwd = 2)
lines(MA20, col = "red", lwd = 2)


# Ajouter une légende
legend("topleft", legend = c("DJI Prices", "10-day MA", "20-day MA"), 
       col = c("black", "blue", "red"), lty = 1, lwd = 2)


# ---- HMM PART ----


prix <- as.numeric(DJI$DJI.Adjusted)  # Pour faciliter les manipulations
Flag <- ifelse(diff(prix) >= 0, "u", "d")  # Création du Flag


# Initialisation des matrices du HMM
A <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
B <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
pi <- c(0.6, 0.4)


DJI_HMM <- HHMM_2State(A, B, pi, Flag, 300)
DJI_Trend <- DJI_HMM$chemin
DJI_Trend <- ifelse(DJI_Trend == 2, -1, 1)  # Transformer en -1 et 1


# Création du graphique des tendances
barplot(DJI_Trend, col = ifelse(DJI_Trend == 1, "green", "red"), border = NA, 
        main = "Market Trend Signal", ylab = "Trend", xlab = "Date", space = 0)


DJI_HMM$A
DJI_HMM$B

