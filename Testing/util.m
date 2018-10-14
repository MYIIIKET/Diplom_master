classdef util
    methods (Static)
        function distances = getDistances(posA, posB)
            A = sum(posA.^2, 2);
            B = sum(posB.^2, 2)';
            AB = posA*posB';
            distances = sqrt(bsxfun(@plus, B, bsxfun(@minus, A, 2*AB)));
        end
        function indeces = getNeighbors(creatureIndex, radius, distances)
            [row, col] = find(distances(creatureIndex,:)>0 & distances(creatureIndex,:)<radius);
            indeces = min([row' col']);
        end
    end
end