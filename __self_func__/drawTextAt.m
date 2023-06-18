function drawTextAt(w,txt,x,y,color)
bRect = Screen('TextBounds',w,txt);
Screen('DrawText',w,txt,x-bRect(3)/2,y-bRect(4)/2,color);
