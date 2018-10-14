function accuracy = Accuracy(f,y,U)
accuracy = sum(sign(f(U))==y(U))/length(U);

end
